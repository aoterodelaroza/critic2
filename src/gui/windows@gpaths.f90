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
  use representations, only: gtrace_point, gtrace_sphere, gtrace_axis, gtrace_axis_both, gtrace_axis_axial,&
     gtrace_axis_plane, gpaths_trace
  use param, only: bohrtoa
  implicit none

  ! starting angles in the plane, in the order of the flux_method_* values
  character(len=*), parameter :: method_names = "Uniform" // c_null_char //&
     "Quotient" // c_null_char // "Dynamical" // c_null_char
  character(len=*), parameter :: method_tooltips = &
     "Directions spread uniformly around the critical point" // c_null_char //&
     "Uniform directions remapped with the quotient of the two Hessian eigenvalues in the plane,&
     & to sample the surface evenly" // c_null_char //&
     "Directions from the linearized dynamical system around the critical point; traces 50 paths&
     & first to estimate the distance to the nearest cage" // c_null_char

  ! a seed CP picked as an atom in the view is the CP at the position
  ! of the atom (the nucleus), within this distance (bohr)
  real*8, parameter :: pick_maxdist = 1d-3

contains

  !> Draw the trace gradient paths window: the form (field, the table
  !> of the seed CPs, the options of each kind of seed CP, style), traced as a blocking job of the main loop
  !> (run_gpaths) into a gradient paths object of the view, whose
  !> editor the Edit Gradient Paths button opens.
  module subroutine draw_gpaths(w)
    use systems, only: sys, sysc, sys_init, ok_system, atlisttype_ncel_frac, atlisttype_nneq
    use representations, only: field_has_cps
    use gui_main, only: ColorGpath, ColorHighlightScene
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_setpos_bottomright,&
       iw_field_combo, iw_calcwidth, iw_combo_simple, iw_inputint, iw_checkbox, iw_coloredit, iw_separatortext,&
       iw_dragfloat_realc
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, ok, ldum, dotrace, doest
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

    ! a new system: the reference field
    if (goodsys) then
       if (.not.w%gp%sys%is(isys) .or. syschanged) then
          call w%gp%sys%set(isys)
          w%gp%fref = field_ref()
          call cancel_pick(w)
       end if
       if (.not.w%gp%fref%ok(sys(isys))) then
          call fref%set(sys(isys),sys(isys)%iref)
          if (.not.fref%same(w%gp%fref)) w%gp%fref = fref
       end if
    end if

    ! a pending pick in the view it was armed on
    if (w%gp%picking > 0) then
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
       call iw_text("Field",alignframe=.true.)
       ifield = w%gp%fref%id
       if (iw_field_combo("##gpfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
          nonestr="<field not available>",sameline=.true.)) then
          call w%gp%fref%set(sys(isys),ifield)
          call cancel_pick(w)
       end if
       call iw_tooltip("Field whose gradient paths are traced",ttshown)
       ifield = w%gp%fref%id
       ok = sys(isys)%goodfield(ifield)
       if (.not.ok) &
          call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       call warn_grid_interp(isys,ifield,"the gradient paths traced on it")
    end if

    if (goodsys .and. ok) then
       ! the seeds: a single shot (the paths up and down from a picked
       ! point), or the selected CPs of the table (the nuclei only, if
       ! there was no search)
       call sync_sel(w,isys)
       call iw_separatortext("Seeds",w%gp%secseeds)
    end if
    if (goodsys .and. ok .and. w%gp%secseeds%open) then
       call iw_text("Single shot",alignframe=.true.)
       call iw_tooltip("Trace the gradient paths up and down from one point",ttshown)
       if (iw_button("Pick##gppickpoint",sameline=.true.,disabled=(w%gp%picking > 0))) &
          call start_pick(w,iview,gppick_point)
       call iw_tooltip("Trace the gradient paths up and down from a point picked in the view",ttshown)
       if (.not.sys(isys)%c%ismolecule) then
          if (iw_checkbox("Symmetry-equivalent paths##gpsym",w%gp%sym,labelleft=.true.)) call sym_changed(w,isys)
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
       if (iw_button("Pick##gppickcp",sameline=.true.,disabled=(w%gp%picking > 0))) &
          call start_pick(w,iview,gppick_cp)
       call iw_tooltip("Pick a critical point in the view and add it to the selection",ttshown)
       call iw_text("(" // string(count(w%gp%sel)) // " selected)",sameline=.true.)
    end if

    ! the paths: the style (more: the options of each kind of seed CP)
    if (goodsys .and. ok) call iw_separatortext("Paths",w%gp%secpaths,hasmore=.true.)
    if (goodsys .and. ok .and. w%gp%secpaths%open) then
       call iw_text("Color",alignframe=.true.)
       ldum = iw_coloredit("##gprgb",rgb=w%gp%rgb,sameline=.true.)
       call iw_tooltip("Color of the new paths",ttshown)
       call iw_text("Radius (Å)",alignframe=.true.,gap=.true.)
       ldum = iw_dragfloat_realc("##gprad",x1=w%gp%rad,speed=0.002_c_float,&
          min=0.005_c_float,max=1._c_float,scale=real(bohrtoa,c_float),decimal=3,&
          sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Radius of the tubes or spheres of the new paths (Å)",ttshown)

       if (w%gp%secpaths%more) then
          call iw_text("NCP",alignframe=.true.)
          call iw_tooltip("The paths traced from the selected nuclei and nuclear critical points,&
             & in directions spread over a sphere",ttshown)
          call sphere_options("n",w%gp%ntheta_n,w%gp%nphi_n)

          call iw_text("BCP",alignframe=.true.)
          call iw_tooltip("The paths traced down from the selected bond critical points, in the plane&
             & perpendicular to the bond path (the interatomic surface; the bond paths are already&
             & traced by the critical point search)",ttshown)
          call iw_text("↓",alignframe=.true.,gap=.true.)
          call plane_options("b",w%gp%npts_b,w%gp%method_b,1_c_int)

          call iw_text("RCP",alignframe=.true.)
          call iw_tooltip("The paths traced from the selected ring critical points: down along the&
             & axis (the ring paths), and up in the perpendicular plane (the ring surface)",ttshown)
          call iw_text("↓",alignframe=.true.,gap=.true.)
          ldum = iw_checkbox("Ring path##gpringpath",w%gp%ringpath,sameline=.true.,labelleft=.true.)
          call iw_tooltip("Trace the two ring paths, down along the axis of the RCP",ttshown)
          call iw_text("↑",alignframe=.true.,gap=.true.)
          call plane_options("r",w%gp%npts_r,w%gp%method_r,0_c_int)

          call iw_text("CCP",alignframe=.true.)
          call iw_tooltip("The paths traced from the selected cage critical points, in directions&
             & spread over a sphere",ttshown)
          call sphere_options("c",w%gp%ntheta_c,w%gp%nphi_c)
       end if
    end if

    if (goodsys .and. ok) then
       ! trace
       formerr = form_error(w,isys)
       dotrace = iw_button("Trace##gptrace",danger=.true.,disabled=(len(formerr) > 0 .or. w%gp%picking > 0))
       call iw_tooltip("Trace the gradient paths of the selected critical points and add them to the&
          & gradient paths object of this field in the view (created if there is none), one group&
          & per critical point (per symmetry-unique one if all the copies are drawn)",ttshown)
       doest = iw_button("Estimate Cost##gpestimate",sameline=.true.,&
          disabled=(len(formerr) > 0 .or. w%gp%picking > 0))
       call iw_tooltip("Estimate how long the trace takes, by tracing the paths of a sample of the&
          & selected critical points",ttshown)
       irep = target_object(w,iview,.false.)
       if (iw_button("Edit Gradient Paths##gpeditobj",sameline=.true.,disabled=(irep == 0))) &
          iaux = stack_create_window(wintype_editrep,.true.,isys=isys,irep=irep,idparent=iview,orraise=-1)
       call iw_tooltip("Open the gradient paths object of this field in the view",ttshown)
       if (dotrace) call request_job(w,iview,gpjob_trace)
       if (doest) call request_job(w,iview,gpjob_estimate)
       if (allocated(w%gp%estimate)) then
          if (.not.same_key(w%gp%estkey,estimate_key(w))) deallocate(w%gp%estimate,w%gp%estkey)
       end if
       if (len(formerr) > 0) then
          call iw_text(formerr,disabled=.true.,wrap=.true.)
       elseif (allocated(w%gp%estimate)) then
          call iw_text(w%gp%estimate,wrap=.true.)
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

  contains
    !> The number of azimuthal and polar angles of a sphere of paths
    !> (ImGui IDs suffixed with id), in the current line.
    subroutine sphere_options(id,ntheta,nphi)
      character(len=*), intent(in) :: id
      integer(c_int), intent(inout) :: ntheta, nphi

      call iw_text("Azimuthal",alignframe=.true.,gap=.true.)
      ldum = iw_inputint("##gpntheta" // id,ntheta,sameline=.true.,width=3)
      call iw_tooltip("Number of azimuthal angles of the starting directions",ttshown)
      call iw_text("Polar",alignframe=.true.,gap=.true.)
      ldum = iw_inputint("##gpnphi" // id,nphi,sameline=.true.,width=3)
      call iw_tooltip("Number of polar angles of the starting directions",ttshown)
      ntheta = max(ntheta,1)
      nphi = max(nphi,0)

    end subroutine sphere_options

    !> The number of paths in the plane (at least nmin) and their
    !> starting angles (ImGui IDs suffixed with id), in the current
    !> line.
    subroutine plane_options(id,npts,method,nmin)
      character(len=*), intent(in) :: id
      integer(c_int), intent(inout) :: npts, method
      integer(c_int), intent(in) :: nmin

      call iw_text("Number",alignframe=.true.,sameline=.true.)
      ldum = iw_inputint("##gpnpts" // id,npts,sameline=.true.,width=3)
      if (nmin == 0) then
         call iw_tooltip("Number of paths in the plane, 0 for none",ttshown)
      else
         call iw_tooltip("Number of paths in the plane",ttshown)
      end if
      npts = max(npts,nmin)
      call iw_text("Angles",alignframe=.true.,gap=.true.)
      call iw_combo_simple("##gpmethod" // id,method_names,method,tooltips=method_tooltips,&
         ttshown=ttshown,sameline=.true.)
      call iw_tooltip("How the starting directions in the plane are chosen",ttshown)

    end subroutine plane_options
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
    if (w%gp%pending_kind == gpjob_estimate) then
       call iw_wait_overlay("Estimating the cost of the trace...",info,.false.)
    else
       call iw_wait_overlay("Tracing gradient paths...",info,.true.)
    end if

  end subroutine block_gpaths

  !> Run the pending trace of the trace gradient paths window, called
  !> by the main loop after the frame with the overlay: trace the
  !> paths of each selected CP (of each symmetry-unique one, if all the
  !> copies are drawn), or the paths up and down from the single-shot
  !> point, and add them, as new groups, to the gradient paths object
  !> of the field in the view the trace was requested in (created if
  !> there is none). The trace can be cancelled with
  !> Esc, and then no paths are added.
  module subroutine run_gpaths(w)
    use systems, only: sys, sys_init
    use gui_main, only: begin_cancellable, end_cancellable
    use global, only: abort_poll, progress_start, progress_step
    class(window), intent(inout), target :: w

    integer :: isys, iview, ifield, irep, i, k, ntraced
    integer, allocatable :: icel(:)
    type(gpaths_trace), allocatable :: t(:)
    type(gpaths_trace) :: tr
    type(rep_gpaths) :: gnew
    character(len=:), allocatable :: errmsg, err1
    integer :: kind
    logical :: cancelled

    ! the job
    kind = w%gp%pending_kind
    w%gp%pending_kind = gpjob_trace

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
    if (kind /= gpjob_point) then
       errmsg = form_error(w,isys)
       if (len(errmsg) > 0) goto 999
    end if

    ! a cost estimate instead of the trace
    if (kind == gpjob_estimate) then
       call estimate_gpaths(w,isys)
       return
    end if

    ! one group per selected CP (in symmetry mode, only the first copy
    ! of each symmetry-unique CP can be selected), except for an RCP:
    ! one per ring path and one for the paths in the plane; traced
    ! apart from the object so a cancel leaves it unchanged; on error,
    ! go on with the next one and report the first error
    ntraced = 0
    errmsg = ""
    gnew%fref = w%gp%fref
    call begin_cancellable()
    if (kind == gpjob_point) then
       ! single shot: one group, the paths up and down from the point;
       ! none leaves it (a CP, e.g. the middle of a symmetric bond) is
       ! an error
       allocate(t(1))
       tr = cp_recipe(w,isys,0)
       tr%x = w%gp%pending_x
       call trace_one(0,0)
       if (ntraced > 0) then
          if (all(t(1)%site(1)%path(1:t(1)%site(1)%npath)%n < 2)) then
             ntraced = 0
             errmsg = "The gradient vanishes at the picked point (a critical point?): no paths were traced"
          end if
       end if
    else
       icel = selected_cps(w,isys)
       allocate(t(3*size(icel,1)))
       call progress_start(size(icel,1),"critical points")
       associate(f => sys(isys)%f(ifield))
         do i = 1, size(icel,1)
            k = icel(i)
            if (abort_poll()) exit

            ! the recipe: the seed and the options of its kind
            tr = cp_recipe(w,isys,f%cpcel(k)%typ)
            tr%x = f%cpcel(k)%x
            if (f%cpcel(k)%typ == 1 .and. tr%iup /= gtrace_axis_plane) then
               call trace_one(gtrace_axis_axial,1)
               call trace_one(gtrace_axis_axial,2)
               if (tr%iup == gtrace_axis_both) call trace_one(gtrace_axis_plane,0)
            else
               call trace_one(tr%iup,0)
            end if
            call progress_step()
         end do
       end associate
    end if
    cancelled = end_cancellable()
    if (cancelled) then
       errmsg = "The trace was cancelled: no paths were added"
       goto 999
    end if

    ! add the groups to the object of the field in the view (created
    ! if there is none), settled for the current field first (an
    ! update could drop stale groups)
    if (ntraced > 0) then
       irep = target_object(w,iview,.true.)
       call win(iview)%sc%rep(irep)%update()
       associate(g => win(iview)%sc%rep(irep)%gpaths)
         g%tfield = gnew%tfield
         g%tcuid = gnew%tcuid
         g%trev = gnew%trev
         g%ttime = gnew%ttime
         do k = 1, ntraced
            call g%add_trace(t(k))
         end do
       end associate
       win(iview)%sc%rep(irep)%shown = .true.
       win(iview)%sc%forcebuildlists = .true.
       call win(iview)%sc%undo_note("Trace gradient paths",.false.)
    end if
    if (len(errmsg) > 0) goto 999
    return

999 continue
    w%errmsg = errmsg

  contains

    !> Trace a new group with recipe tr, with iup and ipath replaced;
    !> on error, keep the first error message and drop the group.
    subroutine trace_one(iup,ipath)
      integer, intent(in) :: iup, ipath

      ntraced = ntraced + 1
      t(ntraced) = tr
      t(ntraced)%iup = iup
      t(ntraced)%ipath = ipath
      call gnew%run_trace(isys,t(ntraced),err1,abortable=.true.)
      if (len(err1) > 0) then
         if (len(errmsg) == 0) errmsg = err1
         ntraced = ntraced - 1
      end if

    end subroutine trace_one

  end subroutine run_gpaths

  !> Stop highlighting the seed CP under the mouse in the table of the
  !> trace gradient paths window.
  module subroutine clear_gpaths_hover(w)
    class(window), intent(inout), target :: w

    call cps_hover_set(w%gp%fref,0,(/0,0/),w%gp%ihover,w%gp%hoverview,w%gp%hoverview_uid)

  end subroutine clear_gpaths_hover

  !xx! private procedures

  !> The table of the CPs of the field, nuclei included, in the colors
  !> of the view iview (the badges): the cell CPs in irow (the
  !> symmetry-unique CPs, as their first cell copy, in symmetry mode).
  !> A click on a row selects or unselects the CP. Arrows in the first
  !> column mark the CPs with paths already known (the bond paths of
  !> the field, and the traced paths of its object in the view), up
  !> (↑) and/or down (↓). Returns in ihover the CP under the mouse:
  !> (symmetry-unique CP, 0) in symmetry mode,
  !> (0, cell CP) otherwise, (0,0) if none.
  subroutine draw_cp_table(w,isys,iview,irow,ihover)
    use systems, only: sys
    use utils, only: iw_text, iw_highlight_selectable, iw_table_column, iw_table_headers_row,&
       iw_calcheight
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    integer, intent(in) :: irow(:)
    integer, intent(out) :: ihover(2)

    integer :: kr, k, iu, icp, j, ncol, nrow
    integer(c_int) :: flags
    logical :: clicked, symm, hov
    logical, allocatable :: known(:,:)
    real*8 :: x(3)
    real(c_float) :: xpos
    character(kind=c_char,len=:), allocatable, target :: str1, suffix
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
      call known_paths(w,isys,iview,symm,known)

      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      str1 = "##tablegpcps" // c_null_char
      sz%x = 0._c_float
      sz%y = iw_calcheight(min(10,nrow+1),0,.false.)
      if (igBeginTable(c_loc(str1),4,flags,sz,0._c_float)) then
         ncol = -1
         call iw_table_column("",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
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
               icp = merge(iu,k,symm)
               hov = .false.
               clicked = .false.
               if (igTableSetColumnIndex(0_c_int)) then
                  ! the row selectable selects the seed
                  xpos = igGetCursorPosX()
                  call iw_text("",alignframe=.true.)
                  hov = iw_highlight_selectable("##gpcprow" // suffix,clicked=clicked,selected=w%gp%sel(k))
                  call igSameLine(0._c_float,0._c_float)
                  call igSetCursorPosX(xpos)

                  ! the paths from this CP that are already known, in
                  ! two fixed slots (up, down; transparent if not known)
                  call arrow("↑",known(1,icp),.false.)
                  call arrow("↓",known(2,icp),.true.)
               end if
               if (igTableSetColumnIndex(1_c_int)) then
                  ! the badge (an atom for the nuclei), which takes the
                  ! mouse from the selectable: it acts as the row
                  call cp_atom_badge(iview,isys,w%gp%fref%id,icp,.not.symm,"gpcpbadge" // suffix)
                  if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                     hov = .true.
                     if (igIsMouseClicked(ImGuiMouseButton_Left,.false._c_bool)) clicked = .true.
                  end if
               end if
               if (hov) ihover = merge((/iu,0/),(/0,k/),symm)

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
               if (igTableSetColumnIndex(2_c_int)) then
                  if (symm) then
                     x = f%cp(iu)%x
                  else
                     x = f%cpcel(k)%x
                  end if
                  if (c%ismolecule) x = (c%x2c(x) + c%molx0) * bohrtoa
                  call iw_text(string(x(1),'f',10,4) // " " // string(x(2),'f',10,4) // " " //&
                     string(x(3),'f',10,4))
               end if
               if (igTableSetColumnIndex(3_c_int)) &
                  call iw_text(string(f%cp(iu)%s%f,'e',decimal=5))
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call igEndTable()
      end if
    end associate

  contains
    !> Draw the arrow str, highlighted if shown, transparent otherwise.
    subroutine arrow(str,shown,sameline)
      character(len=*), intent(in) :: str
      logical, intent(in) :: shown, sameline

      if (shown) then
         call iw_text(str,highlight=.true.,sameline=sameline,nocapture=.true.)
      else
         call iw_text(str,rgba=(/0._c_float,0._c_float,0._c_float,0._c_float/),sameline=sameline,&
            nocapture=.true.)
      end if

    end subroutine arrow

  end subroutine draw_cp_table

  !> The CPs of the field of window w (system isys) with paths
  !> already known: the bond paths of the field (up from a BCP), and
  !> the paths traced in the object of the field in view iview.
  !> known(1,i) if up, known(2,i) if down, for the symmetry-unique CP
  !> i if symm, the cell CP i otherwise.
  subroutine known_paths(w,isys,iview,symm,known)
    use systems, only: sys
    use representations, only: field_has_cps
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(in) :: symm
    logical, allocatable, intent(out) :: known(:,:)

    integer :: irep, it, is, k, iup
    logical :: ud(2)

    associate(f => sys(isys)%f(w%gp%fref%id))
      allocate(known(2,merge(f%ncp,f%ncpcel,symm)))
      known = .false.

      ! the bond paths of the field, traced up if the nuclei are
      ! maxima (as in autocp and write_paths)
      if (field_has_cps(isys,w%gp%fref%id,withpaths=.true.)) then
         iup = merge(1,2,f%typnuc < 0)
         do k = 1, size(known,2)
            known(iup,k) = any(f%cpgp(:,merge(k,f%cpcel(k)%idx,symm))%n > 1)
         end do
      end if

      ! the traced paths, by site
      irep = target_object(w,iview,.false.)
      if (irep == 0) return
      associate(g => win(iview)%sc%rep(irep)%gpaths)
        if (g%traces_stale(isys)) return
        do it = 1, g%ntrace
           if (.not.g%trace(it)%traced) cycle
           do is = 1, g%trace(it)%nsite
              associate(s => g%trace(it)%site(is))
                if (s%icp <= 0 .or. s%npath == 0) cycle
                ud = (/any(s%path(1:s%npath)%iup > 0),any(s%path(1:s%npath)%iup < 0)/)
                if (symm) then
                   known(:,s%icp) = known(:,s%icp) .or. ud
                else
                   do k = 1, s%ncopy
                      known(:,s%icel(k)) = known(:,s%icel(k)) .or. ud
                   end do
                end if
              end associate
           end do
        end do
      end associate
    end associate

  end subroutine known_paths

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

  !> The recipe of the form (system isys) for a seed CP of type typ
  !> (0 = a point, traced up and down), without its position: the
  !> options of its kind, and the style shared by all the groups. An
  !> RCP with neither the ring paths nor paths in the plane is not
  !> traced (form_error).
  function cp_recipe(w,isys,typ) result(t)
    type(window), intent(in), target :: w
    integer, intent(in) :: isys, typ
    type(gpaths_trace) :: t

    select case (typ)
    case (0)
       t%kind = gtrace_point
       t%iup = 0
    case (-3)
       t%kind = gtrace_sphere
       t%ntheta = w%gp%ntheta_n
       t%nphi = w%gp%nphi_n
    case (-1)
       t%kind = gtrace_axis
       t%iup = gtrace_axis_plane
       t%npts = w%gp%npts_b
       t%method = w%gp%method_b
    case (1)
       t%kind = gtrace_axis
       if (w%gp%npts_r == 0) then
          t%iup = gtrace_axis_axial
       elseif (w%gp%ringpath) then
          t%iup = gtrace_axis_both
       else
          t%iup = gtrace_axis_plane
       end if
       t%npts = w%gp%npts_r
       t%method = w%gp%method_r
    case default
       t%kind = gtrace_sphere
       t%ntheta = w%gp%ntheta_c
       t%nphi = w%gp%nphi_c
    end select
    t%sym = sym_mode(w,isys)
    t%rgb = w%gp%rgb
    t%rad = w%gp%rad

  end function cp_recipe

  !> The selected seed cell CPs of the form (system isys): the
  !> selected ones, except those of type 0. The selection must be
  !> valid (form_error).
  function selected_cps(w,isys) result(icel)
    use systems, only: sys
    type(window), intent(in), target :: w
    integer, intent(in) :: isys
    integer, allocatable :: icel(:)

    integer :: k

    associate(f => sys(isys)%f(w%gp%fref%id))
      icel = pack([(k, k=1,f%ncpcel)],w%gp%sel .and. f%cpcel(1:f%ncpcel)%typ /= 0)
    end associate

  end function selected_cps

  !> Estimate the wall time of the trace of the form on system isys:
  !> build the seeds of a sample of the selected CPs and trace them,
  !> one CP at a time as the trace does, until all are done or the
  !> time budget runs out, and scale the elapsed time to all the
  !> selected CPs. The sample follows a golden-ratio sequence, so any
  !> part of it spreads over the whole selection. The first error, if
  !> any, goes to w%errmsg. Sets w%gp%estimate and w%gp%estkey.
  subroutine estimate_gpaths(w,isys)
    use systems, only: sys
    use flux, only: flux_seed, flux_path, flux_trace
    use global, only: prunedist
    use tools_io, only: string
    use utils, only: duration_string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    real*8, parameter :: budget = 2d0 ! seconds
    real*8, parameter :: golden = 0.6180339887498949d0

    type(gpaths_trace) :: t
    type(flux_seed), allocatable :: seeds(:)
    type(flux_path), allocatable :: paths(:)
    integer, allocatable :: icel(:)
    integer :: k, i, ifield, ncp, ndone, npath
    integer*8 :: c0, c1, crate
    real*8 :: tel
    character(len=:), allocatable :: err1

    ifield = w%gp%fref%id
    associate(f => sys(isys)%f(ifield))
      ! the selected seed CPs
      icel = selected_cps(w,isys)
      ncp = size(icel,1)

      ! trace the sample, until the budget runs out
      ndone = 0
      npath = 0
      call system_clock(c0,crate)
      do k = 1, ncp
         i = icel(min(1 + int(ncp * modulo(k * golden,1d0)),ncp))
         t = cp_recipe(w,isys,f%cpcel(i)%typ)
         if (allocated(seeds)) deallocate(seeds)
         call t%make_seeds(isys,ifield,f%cp(f%cpcel(i)%idx)%x,f%cpcel(i)%typ,seeds,err1)
         if (len(err1) > 0 .and. len_trim(w%errmsg) == 0) w%errmsg = err1
         call flux_trace(f,seeds,paths,prunedist)
         ndone = ndone + 1
         npath = npath + size(seeds,1)
         call system_clock(c1)
         if (real(c1-c0,8) / real(crate,8) > budget) exit
      end do
      call system_clock(c1)
      tel = real(c1-c0,8) / real(crate,8) * real(ncp,8) / real(max(ndone,1),8)
    end associate

    w%gp%estimate = "Estimated cost: ~" // duration_string(tel) // " (timed " //&
       string(ndone) // " of " // string(ncp) // " critical points, " // string(npath) // " paths)"
    w%gp%estkey = estimate_key(w)

  end subroutine estimate_gpaths

  !> The inputs of the cost estimate of the form: the options of the
  !> trace, the field, and the selection and when it was made.
  function estimate_key(w) result(key)
    type(window), intent(in), target :: w
    integer*8, allocatable :: key(:)

    key = [int(w%gp%ntheta_n,8), int(w%gp%nphi_n,8), int(w%gp%npts_b,8), int(w%gp%method_b,8),&
       merge(1_8,0_8,w%gp%ringpath), int(w%gp%npts_r,8), int(w%gp%method_r,8), int(w%gp%ntheta_c,8),&
       int(w%gp%nphi_c,8), w%gp%fref%uid, transfer(w%gp%seltime,1_8)]
    if (allocated(w%gp%sel)) key = [key, merge(1_8,0_8,w%gp%sel)]

  end function estimate_key

  !> Whether the two keys are the same.
  function same_key(k1,k2) result(same)
    integer*8, intent(in) :: k1(:), k2(:)
    logical :: same

    same = (size(k1,1) == size(k2,1))
    if (same) same = all(k1 == k2)

  end function same_key

  !> Why the form cannot be traced (empty if it can), for system isys:
  !> no CPs selected, the selection is not for the current field and
  !> CP list, or an RCP is selected with neither the ring paths nor
  !> paths in the plane.
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
    if (.not.w%gp%ringpath .and. w%gp%npts_r == 0) then
       associate(f => sys(isys)%f(w%gp%fref%id))
         if (any(w%gp%sel .and. f%cpcel(1:f%ncpcel)%typ == 1)) err = "Check Ring path or set the RCP Number"
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

  !> Start a pick in view iview of the target kind (gppick_*): a seed
  !> CP (a CP or an atom), or the point to trace the paths from (a
  !> bond, its midpoint, or a point).
  subroutine start_pick(w,iview,kind)
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview, kind

    if (iview < 1 .or. iview > nwin) return
    w%gp%picking = kind
    w%gp%pickview = iview
    call w%gp%pick%arm()
    if (kind == gppick_point) then
       call win(iview)%viewmode_set_forced(vm_pick_bond,"Pick a bond or point to trace the gradient&
          & paths from",w%id,acceptempty=.true.)
    else
       call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick a critical point to trace the gradient&
          & paths from",w%id,acceptempty=.true.)
    end if

  end subroutine start_pick

  !> Forget the pending pick of window w, releasing the view.
  subroutine cancel_pick(w)
    type(window), intent(inout), target :: w

    if (w%gp%picking == 0) return
    if (w%gp%pickview >= 1 .and. w%gp%pickview <= nwin) &
       call win(w%gp%pickview)%viewmode_release_forced(w%id)
    w%gp%picking = 0

  end subroutine cancel_pick

  !> Request a blocking job of kind kind (gpjob_*) for view iview; x
  !> is the single-shot point (crystallographic). The main loop runs
  !> it after the frame with its overlay; a second request in the same
  !> frame is ignored.
  subroutine request_job(w,iview,kind,x)
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview, kind
    real*8, intent(in), optional :: x(3)

    if (.not.w%request_block()) return
    w%gp%pending_kind = kind
    w%gp%pending_view = iview
    if (present(x)) w%gp%pending_x = x
    w%errmsg = ""

  end subroutine request_job

  !> Handle the pending pick in view iview (system isys): the clicked
  !> CP of the field, or the CP at the clicked atom (its nucleus), is
  !> added to the selection; or the single-shot trace from the picked
  !> point is requested.
  subroutine poll_pick(w,isys,iview)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview

    integer :: istat, cpidx(5), k
    real*8 :: xc(3), d

    call view_pick_result(iview,w%id,isys,w%gp%pick,istat,xc,cpidx)
    if (istat == ipick_point .and. w%gp%picking == gppick_point) then
       call request_job(w,iview,gpjob_point,sys(isys)%c%c2x(xc))
    elseif (istat == ipick_point .and. sys(isys)%goodfield(w%gp%fref%id)) then
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
    if (istat /= ipick_pending) w%gp%picking = 0

  end subroutine poll_pick

end submodule gpaths
