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

! Routines for the critical points window (Tools > Critical Points),
! the GUI interface to the AUTO keyword.
submodule (windows) cp
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! names of the seed kinds for the combo, and their AUTO keywords, in
  ! the order of the cpseed_* values
  character(len=*), parameter :: seedkind_names = "Wigner-Seitz (WS)" // c_null_char //&
     "Atom pairs (PAIR)" // c_null_char // "Atom triplets (TRIPLET)" // c_null_char //&
     "Line (LINE)" // c_null_char // "Sphere (SPHERE)" // c_null_char //&
     "Octahedron (OH)" // c_null_char // "Point (POINT)" // c_null_char //&
     "Molecular mesh (MESH)" // c_null_char
  character(len=7), parameter :: seedkind_kw(0:7) = (/"ws     ","pair   ","triplet",&
     "line   ","sphere ","oh     ","point  ","mesh   "/)

  ! maximum WS/OH subdivision level AUTO accepts
  integer, parameter :: maxdepth = 7

contains

  !> Draw the critical points window: the AUTO form for a field of the
  !> system in the anchor view (the Search tab), run as a blocking job
  !> of the main loop (run_cp_pending).
  module subroutine draw_cp(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_begintabitem,&
       iw_setpos_bottomright
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged
    integer :: isys, iview
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: str1

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)

    ! initialize the form on the first pass
    if (w%firstpass) then
       w%cp = cp_state()
       w%errmsg = ""
    end if

    ! a new system: the reference field and AUTO's default seeds for it
    if (goodsys) then
       if (w%cp%isys /= isys) then
          w%cp%isys = isys
          w%cp%ifield = sys(isys)%iref
          call reset_seeds(w,isys)
          if (allocated(w%cp%summary)) deallocate(w%cp%summary)
       end if
    end if

    ! the system
    if (goodsys) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       str1 = "##drawcp_tabbar" // c_null_char
       flags = ImGuiTabBarFlags_None
       if (igBeginTabBar(c_loc(str1),flags)) then
          if (iw_begintabitem("Search##drawcp_searchtab")) then
             call draw_search_tab(w,isys,iview,ttshown)
             call igEndTabItem()
          end if
          call iw_tooltip("Find the critical points of a field (AUTO)",ttshown)
          call igEndTabBar()
       end if
    end if

    ! the error message from the last operation
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    ! close button
    call iw_setpos_bottomright(5,1)
    if (iw_button("Close")) doquit = .true.
    if (iw_close_event(w%focused())) doquit = .true.

    if (doquit) call w%end()

  end subroutine draw_cp

  !> Draw the overlay shown while the blocking job of the critical
  !> points window runs.
  module subroutine block_cp(w)
    use gui_main, only: io, g
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_blank_background
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    integer(c_int) :: flags
    type(ImVec2) :: pos, sz, pivot
    character(kind=c_char,len=:), allocatable, target :: str1, text
    logical(c_bool) :: ldum
    integer :: isys

    call iw_blank_background()

    ! the text of the overlay: what is being done
    text = "...Searching for critical points..."
    isys = w%cp%isys
    if (ok_system(isys,sys_init)) then
       text = text // newline // "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (sys(isys)%goodfield(w%cp%ifield)) &
          text = text // newline // "Field:  " // string(w%cp%ifield) // ": " //&
          trim(sys(isys)%f(w%cp%ifield)%name)
    end if
    if (allocated(w%cp%pending_line)) text = text // newline // "Input:  AUTO " // w%cp%pending_line
    text = text // c_null_char

    ! a centered window, sized for the text: it is only drawn in one
    ! frame, and an auto-resized window is invisible in its first
    pos%x = io%DisplaySize%x * 0.5_c_float
    pos%y = io%DisplaySize%y * 0.5_c_float
    pivot%x = 0.5_c_float
    pivot%y = 0.5_c_float
    call igSetNextWindowPos(pos,0,pivot)
    call igCalcTextSize(sz,c_loc(text),c_null_ptr,.false._c_bool,-1._c_float)
    sz%x = sz%x + 2 * g%Style%WindowPadding%x
    sz%y = sz%y + 2 * g%Style%WindowPadding%y
    call igSetNextWindowSize(sz,0)
    flags = ImGuiWindowFlags_NoDecoration
    flags = ior(flags,ImGuiWindowFlags_NoDocking)
    flags = ior(flags,ImGuiWindowFlags_NoSavedSettings)
    flags = ior(flags,ImGuiWindowFlags_NoFocusOnAppearing)
    flags = ior(flags,ImGuiWindowFlags_NoNav)
    ldum = .true.
    str1 = "##popupwaitcp" // c_null_char
    call igSetNextWindowFocus()
    if (igBegin(c_loc(str1), ldum, flags)) &
       call igTextUnformatted(c_loc(text),c_null_ptr)
    call igEnd()

  end subroutine block_cp

  !> Run the blocking job of the critical points window, called by the
  !> main loop after the frame with the overlay: AUTO on the chosen
  !> field, then the checkpoint, the summary, and the critical points
  !> and gradient paths objects in the view. If AUTO rejects the
  !> options, the field is left as it was.
  module subroutine run_cp_pending(w)
    use autocp, only: autocritic
    use systemmod, only: sy
    use systems, only: sys, sysc, sys_init, ok_system, launch_initialization_thread,&
       kill_initialization_thread, are_threads_running, lastchange_cplist
    use global, only: iunit, iunit_bohr, cp_hdegen
    use tools_io, only: uout
    class(window), intent(inout), target :: w

    integer :: isys, ifield, iref0, iunit0, iview
    real*8 :: hdegen0
    logical :: reinit, ldum, ok
    character(len=:), allocatable :: cpfile, errmsg

    isys = w%cp%isys
    ifield = w%cp%ifield
    iview = w%cp%pending_view
    ok = ok_system(isys,sys_init)
    if (ok) ok = sys(isys)%goodfield(ifield)
    if (ok) ok = allocated(w%cp%pending_line)
    if (.not.ok) then
       w%errmsg = "The system or the field is no longer available"
       return
    end if

    ! stop the initialization threads and connect the system
    reinit = are_threads_running()
    if (reinit) call kill_initialization_thread()
    sy => sys(isys)

    ! AUTO works on the reference field: make the chosen field the
    ! reference for the run (plainly: set_reference would drop the
    ! pointprop lists), with the units the options were written in
    ! (bohr if the current units have no length factor), and restore
    ! the global EPSDEGEN it resets
    iref0 = sys(isys)%iref
    iunit0 = iunit
    hdegen0 = cp_hdegen
    sys(isys)%iref = ifield
    if (iunit /= 1 .and. iunit /= 2) iunit = iunit_bohr
    call autocritic(w%cp%pending_line,ok,clear=w%cp%discard_existing)
    iunit = iunit0
    cp_hdegen = hdegen0
    sys(isys)%iref = iref0

    ! the checkpoint, unless disabled (not AUTO's CHK, which reads an
    ! existing checkpoint instead of searching)
    if (ok .and. .not.w%cp%nochk) then
       cpfile = sys(isys)%f(ifield)%chk_cps_file()
       if (len(cpfile) > 0) then
          call sys(isys)%f(ifield)%write_chk_cps(cpfile,errmsg)
          if (len_trim(errmsg) > 0) then
             write (uout,'("!! Warning !! Could not write the checkpoint ",A,": ",A)') cpfile, errmsg
          else
             write (uout,'("* Critical points written to the checkpoint: ",A/)') cpfile
          end if
       end if
    end if

    ! the output, the threads, and the event
    ldum = read_output_uout(.true.)
    if (reinit) call launch_initialization_thread()
    if (.not.ok) then
       w%errmsg = "AUTO did not run: the options were rejected (see the output console)"
       return
    end if
    call sysc(isys)%post_event(lastchange_cplist)

    ! the summary and the objects in the view
    w%cp%summary = cp_summary(isys,ifield)
    if (iview >= 1 .and. iview <= nwin) then
       if (win(iview)%isopen .and. associated(win(iview)%sc)) then
          if (win(iview)%isys == isys) call win(iview)%sc%show_cps(ifield)
       end if
    end if

  end subroutine run_cp_pending

  !--- private procedures ---

  !> The Search tab: field, seeds, advanced options, and the Run button.
  subroutine draw_search_tab(w,isys,iview,ttshown)
    use gui_main, only: pending_block_window
    use systems, only: sys
    use global, only: iunit
    use utils, only: iw_text, iw_button, iw_tooltip, iw_field_combo, iw_calcwidth,&
       iw_checkbox
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(inout) :: ttshown

    integer :: i, idel
    integer(c_int) :: ifield
    logical :: ldum
    character(len=:), allocatable :: xunit, errrun

    if (sys(isys)%c%ismolecule) then
       xunit = " (Å)"
    else
       xunit = " (fractional)"
    end if

    ! field
    call iw_text("Field",highlight=.true.)
    call igSameLine(0._c_float,-1._c_float)
    ifield = w%cp%ifield
    if (iw_field_combo("##cpfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
       nonestr="<field not available>")) w%cp%ifield = ifield
    call iw_tooltip("Field whose critical points are searched",ttshown)
    if (.not.sys(isys)%goodfield(w%cp%ifield)) then
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       return
    end if

    ! the seeds
    call iw_text("Seeds",highlight=.true.)
    call iw_tooltip("Starting points of the critical point searches. Each entry generates a set&
       & of points; the default is the seeding AUTO uses for this system",ttshown)
    idel = 0
    do i = 1, size(w%cp%seed)
       if (draw_seed(w,isys,i,xunit,ttshown)) idel = i
    end do
    if (idel > 0) w%cp%seed = [w%cp%seed(:idel-1), w%cp%seed(idel+1:)]
    if (iw_button("Add seed##cpseedadd")) call add_seed(w,isys,cpseed_point)
    call iw_tooltip("Add a seed (a point at the center of the system; change its kind above)",ttshown)
    if (iw_button("Reset seeds##cpseedreset",sameline=.true.)) call reset_seeds(w,isys)
    call iw_tooltip("Go back to the seeding AUTO uses by default for this system",ttshown)

    ! advanced options
    call draw_advanced(w,isys,xunit,ttshown)

    ! run options
    call iw_text("Run",highlight=.true.)
    ldum = iw_checkbox("Discard the existing critical points##cpdiscardexist",w%cp%discard_existing)
    call iw_tooltip("Start from the nuclei only. Otherwise, the new critical points are added&
       & to those the field already has (as the AUTO keyword does)",ttshown)
    ldum = iw_checkbox("Do not write the checkpoint##cpnochk",w%cp%nochk)
    call iw_tooltip("Do not save the critical points to the checkpoint file (<field file>.chk_cps),&
       & which the GUI reads the next time the field is loaded",ttshown)

    ! the form must make valid AUTO options
    errrun = form_error(w)
    if (iw_button("Run##cprun",disabled=(len(errrun) > 0))) then
       w%cp%pending_line = auto_options(w,isys,iunit)
       w%cp%pending_view = iview
       w%errmsg = ""
       pending_block_window = w%id
    end if
    call iw_tooltip("Search for the critical points (the calculation blocks the interface&
       & until it finishes; the output goes to the output console)",ttshown)
    if (len(errrun) > 0) call iw_text(errrun,danger=.true.,sameline=.true.)

    ! the result of the last run
    if (allocated(w%cp%summary)) call iw_text(w%cp%summary,wrap=.true.)

  end subroutine draw_search_tab

  !> Draw the widgets of seed i of window w: its kind, a button to
  !> remove it (returns true if pressed), and the parameters that kind
  !> uses. Choosing a kind gives its parameters usable values.
  function draw_seed(w,isys,i,xunit,ttshown) result(del)
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_inputint, iw_dragfloat_real8,&
       iw_close_button
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, i
    character(len=*), intent(in) :: xunit
    logical, intent(inout) :: ttshown
    logical :: del

    character(len=:), allocatable :: suf
    logical :: ldum, ch

    suf = "##cpseed" // string(i)
    call iw_text(string(i) // ".",alignframe=.true.)
    call iw_combo_simple("##cpseedkind" // suf,seedkind_names,w%cp%seed(i)%typ,sameline=.true.,&
       changed=ch)
    call iw_tooltip("Kind of seeding",ttshown)
    if (ch) call kind_defaults(w,isys,i)
    call igSameLine(0._c_float,-1._c_float)
    del = iw_close_button("##cpseeddel" // suf)
    call iw_tooltip("Remove this seed",ttshown)

    call igIndent(0._c_float)
    associate(s => w%cp%seed(i))
      select case(s%typ)
      case (cpseed_ws,cpseed_oh)
         ldum = iw_inputint("Subdivision level" // suf // "dep",s%depth,width=6)
         s%depth = max(0,min(s%depth,maxdepth))
         call iw_tooltip("Number of recursive subdivisions of the region (0 to 7; more seeds&
            & at each level)",ttshown)
         call x0_widget("Center" // xunit)
         if (s%typ == cpseed_ws) then
            ldum = iw_dragfloat_real8("Radius (Å, 0 = whole cell)" // suf // "rad",x1=s%rad,&
               speed=0.01d0,min=0d0,max=100d0,decimal=3)
            call iw_tooltip("Scale the Wigner-Seitz cell to this radius around the center&
               & (0 = the whole cell)",ttshown)
         else
            ldum = iw_dragfloat_real8("Radius (Å)" // suf // "rad",x1=s%rad,speed=0.01d0,&
               min=0.01d0,max=100d0,decimal=3)
            ldum = iw_inputint("Radial points" // suf // "nr",s%nr,width=6)
            s%nr = max(s%nr,1)
         end if
      case (cpseed_sphere)
         call x0_widget("Center" // xunit)
         ldum = iw_dragfloat_real8("Radius (Å)" // suf // "rad",x1=s%rad,speed=0.01d0,&
            min=0.01d0,max=100d0,decimal=3)
         ldum = iw_inputint("Polar points" // suf // "nth",s%ntheta,width=6)
         ldum = iw_inputint("Azimuthal points" // suf // "nph",s%nphi,width=6)
         ldum = iw_inputint("Radial points" // suf // "nr",s%nr,width=6)
         s%ntheta = max(s%ntheta,1)
         s%nphi = max(s%nphi,1)
         s%nr = max(s%nr,1)
      case (cpseed_pair,cpseed_triplet)
         ldum = iw_dragfloat_real8("Maximum distance (Å)" // suf // "dist",x1=s%dist,speed=0.01d0,&
            min=0d0,max=100d0,decimal=3)
         call iw_tooltip("Only atoms closer than this distance are combined",ttshown)
         if (s%typ == cpseed_pair) then
            ldum = iw_inputint("Points per pair" // suf // "npts",s%npts,width=6)
            s%npts = max(s%npts,1)
            call iw_tooltip("Number of seeds on the segment between the two atoms",ttshown)
         end if
      case (cpseed_line)
         call x0_widget("Start" // xunit)
         ldum = iw_dragfloat_real8("End" // xunit // suf // "x1",x3=s%x1,speed=0.001d0,decimal=4)
         ldum = iw_inputint("Points" // suf // "npts",s%npts,width=6)
         s%npts = max(s%npts,2)
      case (cpseed_point)
         call x0_widget("Position" // xunit)
      case (cpseed_mesh)
         call iw_text("The points of the molecular integration mesh",disabled=.true.)
      end select
    end associate
    call igUnindent(0._c_float)

  contains
    subroutine x0_widget(label)
      character(len=*), intent(in) :: label
      ldum = iw_dragfloat_real8(label // suf // "x0",x3=w%cp%seed(i)%x0,speed=0.001d0,decimal=4)
    end subroutine x0_widget
  end function draw_seed

  !> The advanced options of AUTO, in a collapsed tree node; each one
  !> is used only if its checkbox is set (otherwise, AUTO's default).
  subroutine draw_advanced(w,isys,xunit,ttshown)
    use utils, only: iw_tooltip, iw_checkbox, iw_inputtext, iw_dragfloat_real8,&
       iw_combo_simple
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    character(len=*), intent(in) :: xunit
    logical, intent(inout) :: ttshown

    character(len=:,kind=c_char), allocatable, target :: strad
    logical :: ldum, ch

    strad = "Advanced options##cpadvanced" // c_null_char
    if (igTreeNodeEx_Str(c_loc(strad),ImGuiTreeNodeFlags_None)) then
       ldum = iw_checkbox("##cpusegradeps",w%cp%use_gradeps)
       ldum = iw_inputtext("Gradient norm (GRADEPS)##cpgradeps",bufsize=31,textf=w%cp%gradeps,&
          width=12,sameline=.true.)
       call iw_tooltip("Maximum gradient norm of a critical point (default 1e-12)",ttshown)

       ldum = iw_checkbox("##cpusecpeps",w%cp%use_cpeps)
       ldum = iw_dragfloat_real8("CP distance (CPEPS, Å)##cpcpeps",x1=w%cp%cpeps,speed=0.001d0,&
          min=0d0,max=10d0,decimal=4,sameline=.true.)
       call iw_tooltip("Two critical points closer than this are the same (default 0.0053 Å)",ttshown)

       ldum = iw_checkbox("##cpusenuceps",w%cp%use_nuceps)
       ldum = iw_dragfloat_real8("Nucleus distance (NUCEPS, Å)##cpnuceps",x1=w%cp%nuceps,speed=0.001d0,&
          min=0d0,max=10d0,decimal=4,sameline=.true.)
       call iw_tooltip("Discard critical points closer than this to a nucleus (default 0.053 Å;&
          & larger for grids)",ttshown)

       ldum = iw_checkbox("##cpusenucepsh",w%cp%use_nucepsh)
       ldum = iw_dragfloat_real8("Hydrogen distance (NUCEPSH, Å)##cpnucepsh",x1=w%cp%nucepsh,&
          speed=0.001d0,min=0d0,max=10d0,decimal=4,sameline=.true.)
       call iw_tooltip("Same, for hydrogen nuclei (default 0.106 Å; larger for grids)",ttshown)

       ldum = iw_checkbox("##cpuseepsdegen",w%cp%use_epsdegen)
       ldum = iw_inputtext("Degenerate eigenvalue (EPSDEGEN)##cpepsdegen",bufsize=31,textf=w%cp%epsdegen,&
          width=12,sameline=.true.)
       call iw_tooltip("Hessian eigenvalues smaller than this make a degenerate critical point&
          & (default 1e-8)",ttshown)

       ldum = iw_inputtext("Discard (DISCARD)##cpdiscard",bufsize=1023,textf=w%cp%discard,width=30)
       call iw_tooltip("Discard the critical points where this expression is not zero, e.g.&
          & $0 < 1e-3 (empty = none)",ttshown)

       call iw_combo_simple("Clip (CLIP)##cpclip","None" // c_null_char // "Box" // c_null_char //&
          "Sphere" // c_null_char,w%cp%iclip,changed=ch)
       call iw_tooltip("Only use the seeds inside this region",ttshown)
       if (ch) call clip_defaults(w,isys)
       if (w%cp%iclip == 1) then
          ldum = iw_dragfloat_real8("Box corner 1" // xunit // "##cpclipx0",x3=w%cp%clipx0,&
             speed=0.001d0,decimal=4)
          ldum = iw_dragfloat_real8("Box corner 2" // xunit // "##cpclipx1",x3=w%cp%clipx1,&
             speed=0.001d0,decimal=4)
       elseif (w%cp%iclip == 2) then
          ldum = iw_dragfloat_real8("Center" // xunit // "##cpclipx0s",x3=w%cp%clipx0,&
             speed=0.001d0,decimal=4)
          ldum = iw_dragfloat_real8("Radius (Å)##cpcliprad",x1=w%cp%cliprad,speed=0.01d0,&
             min=0d0,max=100d0,decimal=3)
       end if
       call igTreePop()
    end if

  end subroutine draw_advanced

  !> Reset the seeds of the window to AUTO's default for system isys:
  !> Wigner-Seitz for crystals, atom pairs for molecules.
  subroutine reset_seeds(w,isys)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    if (allocated(w%cp%seed)) deallocate(w%cp%seed)
    allocate(w%cp%seed(0))
    if (sys(isys)%c%ismolecule) then
       call add_seed(w,isys,cpseed_pair)
    else
       call add_seed(w,isys,cpseed_ws)
    end if

  end subroutine reset_seeds

  !> Append a seed of kind typ with the defaults of that kind.
  subroutine add_seed(w,isys,typ)
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    integer, intent(in) :: typ

    w%cp%seed = [w%cp%seed, cp_seed_ui(typ=typ)]
    call kind_defaults(w,isys,size(w%cp%seed))

  end subroutine add_seed

  !> Give seed i the positions a seed of its kind starts with: the
  !> center of the system (the cell origin for a Wigner-Seitz seed in
  !> a crystal, as AUTO) and, for the kinds that need them, a radius
  !> and a line end.
  subroutine kind_defaults(w,isys,i)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, i

    logical :: ismol

    ismol = sys(isys)%c%ismolecule
    associate(s => w%cp%seed(i))
      s%x0 = system_center(isys)
      if (s%typ == cpseed_ws .and. .not.ismol) s%x0 = 0d0
      if (s%typ == cpseed_line) then
         s%x1 = s%x0
         if (ismol) then
            s%x1(1) = s%x1(1) + 1d0
         else
            s%x1(1) = s%x1(1) + 0.5d0
         end if
      end if
      if ((s%typ == cpseed_sphere .or. s%typ == cpseed_oh) .and. s%rad <= 0d0) s%rad = 1d0
    end associate

  end subroutine kind_defaults

  !> Give the clip region of window w starting values for system isys:
  !> the whole cell or a sphere around its center (crystals), or a box
  !> or sphere of 5 Å around the center of the molecule.
  subroutine clip_defaults(w,isys)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    real*8 :: xc(3)

    xc = system_center(isys)
    w%cp%clipx0 = xc
    w%cp%cliprad = 5d0
    if (w%cp%iclip == 1) then
       if (sys(isys)%c%ismolecule) then
          w%cp%clipx0 = xc - 5d0
          w%cp%clipx1 = xc + 5d0
       else
          w%cp%clipx0 = 0d0
          w%cp%clipx1 = 1d0
       end if
    end if

  end subroutine clip_defaults

  !> The center of system isys, in the coordinates of the form: the
  !> center of the cell (fractional) for crystals, the centroid of the
  !> atoms in the input frame (Å) for molecules.
  function system_center(isys) result(xc)
    use systems, only: sys
    integer, intent(in) :: isys
    real*8 :: xc(3)

    integer :: i

    associate(c => sys(isys)%c)
      if (c%ismolecule .and. c%ncel > 0) then
         xc = 0d0
         do i = 1, c%ncel
            xc = xc + c%atcel(i)%r
         end do
         xc = (xc / c%ncel + c%molx0) * bohrtoa
      else
         xc = 0.5d0
      end if
    end associate

  end function system_center

  !> Why the form of window w would not make valid AUTO options, or an
  !> empty string if it does.
  function form_error(w) result(errmsg)
    type(window), intent(in) :: w
    character(len=:), allocatable :: errmsg

    errmsg = ""
    if (size(w%cp%seed) == 0) then
       errmsg = "Add a seed"
    elseif (w%cp%use_gradeps .and. .not.isnumber(w%cp%gradeps)) then
       errmsg = "GRADEPS is not a number"
    elseif (w%cp%use_epsdegen .and. .not.isnumber(w%cp%epsdegen)) then
       errmsg = "EPSDEGEN is not a number"
    elseif (index(w%cp%discard,'"') > 0) then
       errmsg = "The DISCARD expression cannot contain double quotes"
    end if

  contains
    !> Whether str is a single real number (a list-directed read
    !> alone would accept trailing words).
    function isnumber(str) result(ok)
      character(len=*), intent(in) :: str
      logical :: ok
      real*8 :: x
      integer :: ios
      ok = (len_trim(str) > 0)
      if (ok) ok = (verify(trim(adjustl(str)),"0123456789.+-eEdD") == 0)
      if (ok) then
         read (str,*,iostat=ios) x
         ok = (ios == 0)
      end if
    end function isnumber
  end function form_error

  !> The AUTO options (the text after the keyword) for the form of
  !> window w and system isys, with the lengths in the input units
  !> iunit (bohr if they have no length factor, as run_cp_pending runs
  !> AUTO): positions fractional for crystals and Cartesian in the
  !> input frame for molecules. CPEPS, NUCEPS, and NUCEPSH are always
  !> in bohr, as AUTO reads them.
  function auto_options(w,isys,iunit) result(line)
    use systems, only: sys
    use global, only: dunit0
    use tools_io, only: string
    type(window), intent(in) :: w
    integer, intent(in) :: isys, iunit
    character(len=:), allocatable :: line

    integer :: i
    logical :: ismol
    real*8 :: fl

    ismol = sys(isys)%c%ismolecule
    fl = 1d0 / bohrtoa ! Å to the input units
    if (iunit == 1 .or. iunit == 2) fl = fl * dunit0(iunit)

    ! seeds
    line = ""
    do i = 1, size(w%cp%seed)
       associate(s => w%cp%seed(i))
         line = line // " seed " // trim(seedkind_kw(s%typ))
         select case(s%typ)
         case (cpseed_ws,cpseed_oh)
            line = line // " depth " // string(s%depth) // xstr(" x0",s%x0)
            if (s%rad > 0d0) line = line // " radius " // rstr(s%rad * fl)
            if (s%typ == cpseed_oh) line = line // " nr " // string(s%nr)
         case (cpseed_sphere)
            line = line // xstr(" x0",s%x0) // " radius " // rstr(s%rad * fl) //&
               " ntheta " // string(s%ntheta) // " nphi " // string(s%nphi) // " nr " // string(s%nr)
         case (cpseed_pair)
            line = line // " dist " // rstr(s%dist * fl) // " npts " // string(s%npts)
         case (cpseed_triplet)
            line = line // " dist " // rstr(s%dist * fl)
         case (cpseed_line)
            line = line // xstr(" x0",s%x0) // xstr(" x1",s%x1) // " npts " // string(s%npts)
         case (cpseed_point)
            line = line // xstr(" x0",s%x0)
         end select
       end associate
    end do

    ! advanced options
    if (w%cp%use_gradeps) line = line // " gradeps " // trim(adjustl(w%cp%gradeps))
    if (w%cp%use_cpeps) line = line // " cpeps " // rstr(w%cp%cpeps / bohrtoa)
    if (w%cp%use_nuceps) line = line // " nuceps " // rstr(w%cp%nuceps / bohrtoa)
    if (w%cp%use_nucepsh) line = line // " nucepsh " // rstr(w%cp%nucepsh / bohrtoa)
    if (w%cp%use_epsdegen) line = line // " epsdegen " // trim(adjustl(w%cp%epsdegen))
    if (len_trim(w%cp%discard) > 0) line = line // ' discard "' // trim(w%cp%discard) // '"'
    if (w%cp%iclip == 1) then
       line = line // " clip cube" // xstr("",w%cp%clipx0) // xstr("",w%cp%clipx1)
    elseif (w%cp%iclip == 2) then
       line = line // " clip sphere" // xstr("",w%cp%clipx0) // " " // rstr(w%cp%cliprad * fl)
    end if
    line = line(2:)

  contains
    !> A real number in fixed notation, 10 decimals without the
    !> trailing zeros.
    function rstr(x) result(str)
      real*8, intent(in) :: x
      character(len=:), allocatable :: str
      integer :: n
      str = string(x,'f',decimal=10)
      n = len(str)
      do while (n > 1 .and. str(n:n) == "0")
         n = n - 1
      end do
      if (str(n:n) == ".") n = n - 1
      str = str(1:n)
    end function rstr
    !> key followed by a position: fractional for crystals, the
    !> input-frame Cartesian position (Å, converted to the input
    !> units) for molecules.
    function xstr(key,x) result(str)
      character(len=*), intent(in) :: key
      real*8, intent(in) :: x(3)
      character(len=:), allocatable :: str
      real*8 :: fx
      fx = 1d0
      if (ismol) fx = fl
      str = key // " " // rstr(x(1)*fx) // " " // rstr(x(2)*fx) // " " // rstr(x(3)*fx)
    end function xstr
  end function auto_options

  !> One-line summary of the critical points of field ifield of system
  !> isys: the number of each type in the cell (or molecule) and the
  !> Morse (crystals) or Poincare-Hopf (molecules) sum.
  function cp_summary(isys,ifield) result(str)
    use systems, only: sys
    use tools_io, only: string
    integer, intent(in) :: isys, ifield
    character(len=:), allocatable :: str

    integer :: i, it, nt(0:3)

    nt = 0
    associate(f => sys(isys)%f(ifield))
      do i = 1, f%ncp
         it = f%cp(i)%typind
         if (it >= 0 .and. it <= 3) nt(it) = nt(it) + f%cp(i)%mult
      end do
    end associate
    str = "Critical points (n|b|r|c): " // string(nt(0)) // " | " // string(nt(1)) // " | " //&
       string(nt(2)) // " | " // string(nt(3))
    if (sys(isys)%c%ismolecule) then
       str = str // "; Poincare-Hopf sum: " // string(nt(0)-nt(1)+nt(2)-nt(3))
    else
       str = str // "; Morse sum: " // string(nt(0)-nt(1)+nt(2)-nt(3))
    end if

  end function cp_summary

end submodule cp
