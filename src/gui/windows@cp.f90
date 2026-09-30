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
  ! description of each seed kind, for the tooltips of the kind combo
  character(len=*), parameter :: seedkind_tooltips = &
     "The irreducible part of the Wigner-Seitz cell, split into tetrahedra and subdivided&
     & recursively (more seeds at each level). AUTO's default for crystals." // c_null_char //&
     "Points on the segment between every pair of atoms closer than the maximum distance.&
     & AUTO's default for molecules; good at finding bond critical points." // c_null_char //&
     "The centroid of every triplet of atoms closer to each other than the maximum distance;&
     & aimed at ring critical points." // c_null_char //&
     "Evenly spaced points on the segment between two positions." // c_null_char //&
     "Points on concentric spherical shells around a center (polar x azimuthal x radial)." //&
     c_null_char //&
     "Points on concentric octahedral shells around a center, the octahedron subdivided&
     & recursively." // c_null_char //&
     "A single starting point." // c_null_char //&
     "The points of the molecular integration mesh: dense, and expensive in large systems." //&
     c_null_char

  ! maximum WS/OH subdivision level AUTO accepts
  integer, parameter :: maxdepth = 7

  ! the global state a run of AUTO or CPREPORT from this window changes
  ! (auto_enter/auto_leave): the reference field of the system, the
  ! input units, and EPSDEGEN
  type auto_context
     integer :: iref = 0
     integer :: iunit = 0
     real*8 :: hdegen = 0d0
  end type auto_context

  ! seed preview: at most this many seeds are drawn, as spheres of
  ! this radius (bohr) and color
  integer, parameter :: maxseedshow = 20000
  real*8, parameter :: seed_rad = 0.08d0 / bohrtoa
  real(c_float), parameter :: seed_rgb(3) = (/0.45_c_float,0.45_c_float,0.45_c_float/)

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

    logical :: doquit, goodsys, syschanged, tabopen
    integer :: isys, iview, ihover(2)
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

    ! a new system: the reference field and AUTO's default seeds for
    ! it, and the results table caches are stale
    if (goodsys) then
       if (w%cp%isys /= isys .or. syschanged) then
          w%cp%isys = isys
          w%cp%ifield = sys(isys)%iref
          call reset_seeds(w,isys)
          if (allocated(w%cp%summary)) deallocate(w%cp%summary)
          ! the default export file: <root>.cps.cif (<root>.cif could be the source file)
          w%okfile = okfile_default(isys,"structure","cps.cif")
          w%cp%tfield = -1
       end if
    end if

    ! a pending pick of the point to add a CP from, in the view it was
    ! armed on; dropped, releasing that view, if the window closes, the
    ! system is not ready, or the window moved to another view
    if (w%cp%picking) then
       if (goodsys .and. .not.doquit .and. iview == w%cp%pickview) then
          call poll_add_pick(w,isys,iview)
       else
          if (w%cp%pickview >= 1 .and. w%cp%pickview <= nwin) &
             call win(w%cp%pickview)%viewmode_release_forced(w%id)
          w%cp%picking = .false.
       end if
    end if

    ! the system
    ihover = 0
    if (goodsys) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       str1 = "##drawcp_tabbar" // c_null_char
       flags = ImGuiTabBarFlags_None
       if (igBeginTabBar(c_loc(str1),flags)) then
          ! each tooltip goes right after its tab: after the tab's
          ! contents, it would apply to their last widget
          tabopen = iw_begintabitem("Search##drawcp_searchtab")
          call iw_tooltip("Find the critical points of a field (AUTO)",ttshown)
          if (tabopen) then
             call draw_search_tab(w,isys,iview,ttshown)
             call igEndTabItem()
          end if
          tabopen = iw_begintabitem("Results##drawcp_resultstab")
          call iw_tooltip("The critical points of the field",ttshown)
          if (tabopen) then
             call draw_results_tab(w,isys,iview,ihover,ttshown)
             call igEndTabItem()
          end if
          tabopen = iw_begintabitem("Export##drawcp_exporttab")
          call iw_tooltip("Write the critical points to a file",ttshown)
          if (tabopen) then
             call draw_export_tab(w,isys,iview,ttshown)
             call igEndTabItem()
          end if
          call igEndTabBar()
       end if
    end if

    ! the seeds of the form, shown in the view
    if (goodsys .and. .not.doquit .and. w%cp%showseeds) call show_seed_preview(w,isys)

    ! the CP under the mouse in the results table is highlighted in the view
    if (goodsys) then
       call set_cp_hover(w,iview,ihover)
    else
       call w%clear_cp_hover()
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
    select case (w%cp%pending_kind)
    case (cpjob_add)
       text = "...Searching for a critical point..."
    case (cpjob_delete)
       text = "...Deleting critical points and tracing the bond paths..."
    case (cpjob_export)
       text = "...Writing the critical points..."
    case default
       text = "...Searching for critical points..."
    end select
    isys = w%cp%isys
    if (ok_system(isys,sys_init)) then
       text = text // newline // "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (sys(isys)%goodfield(w%cp%ifield)) &
          text = text // newline // "Field:  " // string(w%cp%ifield) // ": " //&
          trim(sys(isys)%f(w%cp%ifield)%name)
    end if
    if (allocated(w%cp%pending_line)) then
       if (w%cp%pending_kind == cpjob_export) then
          text = text // newline // "Input:  CPREPORT " // w%cp%pending_line
       elseif (w%cp%pending_kind /= cpjob_delete) then
          text = text // newline // "Input:  AUTO " // w%cp%pending_line
       end if
    end if
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
  !> main loop after the frame with the overlay, on the chosen field: a
  !> search (AUTO with the options of the Search tab), the addition of
  !> a CP (AUTO from one point, appending), or the deletion of the
  !> selected CPs (and the bond graph traced again). Then the
  !> checkpoint, the summary, and the critical points and gradient
  !> paths objects in the view. If AUTO rejects the options, the field
  !> is left as it was.
  module subroutine run_cp_pending(w)
    use autocp, only: autocritic, autocritic_graph, cpreport
    use systems, only: sys, sysc, sys_init, ok_system, launch_initialization_thread,&
       kill_initialization_thread, are_threads_running, lastchange_cplist
    use tools_io, only: uout, string
    class(window), intent(inout), target :: w

    integer :: isys, ifield, iview, ncp0, kind
    type(auto_context) :: ctx
    logical :: reinit, ldum, ok
    character(len=:), allocatable :: cpfile, errmsg

    isys = w%cp%isys
    ifield = w%cp%ifield
    iview = w%cp%pending_view
    kind = w%cp%pending_kind
    ok = ok_system(isys,sys_init)
    if (ok) ok = sys(isys)%goodfield(ifield)
    if (ok) then
       if (kind == cpjob_delete) then
          ok = allocated(w%cp%sel)
          if (ok) ok = (size(w%cp%sel) == sys(isys)%f(ifield)%ncp)
       else
          ok = allocated(w%cp%pending_line)
       end if
    end if
    if (.not.ok) then
       w%errmsg = "The system or the field is no longer available"
       return
    end if

    ! stop the initialization threads, and run on the chosen field
    reinit = are_threads_running()
    if (reinit) call kill_initialization_thread()
    call auto_enter(isys,ifield,ctx)
    ncp0 = sys(isys)%f(ifield)%ncp
    if (kind == cpjob_delete) then
       write (uout,'("* Deleting ",A," critical points (GUI) and tracing the bond paths again")') &
          string(count(w%cp%sel))
       call sys(isys)%f(ifield)%delete_cps(w%cp%sel)
       call autocritic_graph()
       write (uout,*)
       ok = .true.
    elseif (kind == cpjob_export) then
       call cpreport(w%cp%pending_line)
       write (uout,'("* Critical points written to: ",A/)') trim(w%okfile)
       ok = .true.
    else
       call autocritic(w%cp%pending_line,ok,clear=(kind == cpjob_search .and. w%cp%discard_existing))
    end if
    call auto_leave(isys,ctx)

    ! the checkpoint, unless disabled (not AUTO's CHK, which reads an
    ! existing checkpoint instead of searching)
    if (ok .and. .not.w%cp%nochk .and. kind /= cpjob_export) then
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

    ! an export leaves the CP list as it was
    if (kind == cpjob_export) then
       call okfile_save_dir(w%okfile)
       return
    end if
    call sysc(isys)%post_event(lastchange_cplist)

    ! a search from one point that found nothing new
    if (kind == cpjob_add .and. sys(isys)%f(ifield)%ncp == ncp0) &
       w%errmsg = "No new critical point: the search from that point found none, or one already in the list"

    ! the summary and the objects in the view
    w%cp%summary = cp_summary(isys,ifield)
    if (iview >= 1 .and. iview <= nwin) then
       if (win(iview)%isopen .and. associated(win(iview)%sc)) then
          if (win(iview)%isys == isys) call win(iview)%sc%show_cps(ifield)
       end if
    end if

  end subroutine run_cp_pending

  !> Stop highlighting in the view the CP under the mouse in the
  !> results table of the critical points window.
  module subroutine clear_cp_hover(w)
    class(window), intent(inout), target :: w

    call set_cp_hover(w,0,(/0,0/))

  end subroutine clear_cp_hover

  !--- private procedures ---

  !> The Results tab: the table of the critical points of the field,
  !> for the symmetry-unique CPs or the cell CPs. Returns in ihover the
  !> CP under the mouse (as rep_cps%ihover).
  subroutine draw_results_tab(w,isys,iview,ihover,ttshown)
    use systems, only: sys, sysc, atlisttype_nneq, atlisttype_ncel_frac
    use representations, only: field_has_cps
    use gui_main, only: g, ColorHighlightScene
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_calcheight, iw_table_column,&
       iw_table_headers_row, iw_highlight_selectable, iw_atom_button
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    integer, intent(inout) :: ihover(2)
    logical, intent(inout) :: ttshown

    integer :: n, k, ifield, ncol, nrow, ihnuc(2)
    integer(c_int) :: itable, flags
    logical :: ch, cell, ismol
    character(kind=c_char,len=:), allocatable, target :: str1
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    type(ImVec2) :: sz

    ifield = w%cp%ifield
    if (.not.sys(isys)%goodfield(ifield)) then
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       return
    end if
    call iw_text("Field",highlight=.true.)
    call iw_text(string(ifield) // ": " // trim(sys(isys)%f(ifield)%name),sameline=.true.)
    if (.not.field_has_cps(isys,ifield)) then
       call iw_text("This field has no critical points other than the nuclei (search for them in&
          & the Search tab, or load a checkpoint that has them)",disabled=.true.,wrap=.true.)
       return
    end if
    call update_table_caches(w,isys)
    ismol = sys(isys)%c%ismolecule
    ihnuc = 0

    associate(f => sys(isys)%f(ifield))
      ! the list: the symmetry-unique CPs, or the cell CPs if there are
      ! more (the combo items are in the order of the tablecell values)
      cell = .false.
      if (f%ncpcel > f%ncp) then
         itable = int(w%cp%tablecell,c_int)
         call iw_combo_simple("Critical point list##cptablecombo","Symmetry-unique" // c_null_char //&
            "Cell" // c_null_char,itable,changed=ch)
         call iw_tooltip("List the symmetry-unique critical points or every critical point in the cell",&
            ttshown)
         if (ch) w%cp%tablecell = int(itable)
         cell = (w%cp%tablecell == 1)
      end if
      n = merge(f%ncpcel,f%ncp,cell)
      call iw_text(w%cp%summary,wrap=.true.)

      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollX)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      str1 = "##tablecpresults" // c_null_char
      sz%x = 0._c_float
      ! rows one frame high (badges) with cell padding, and the
      ! horizontal scrollbar
      nrow = min(16,n+1)
      sz%y = iw_calcheight(nrow,0,.false.) + nrow * 2 * g%Style%CellPadding%y +&
         g%Style%ScrollbarSize
      if (igBeginTable(c_loc(str1),10,flags,sz,0._c_float)) then
         ncol = -1
         call iw_table_column("CP",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         if (ismol) then
            call iw_table_column("Position (Å)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         else
            call iw_table_column("Position (fractional)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         end if
         call iw_table_column("Site",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Field",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("|Gradient|",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Laplacian",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("End 1",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("End 2",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Path (Å)",icol=ncol,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Ellipticity",icol=ncol,flags=ImGuiTableColumnFlags_WidthStretch)
         call iw_table_headers_row(freezetop=.true.,autofit=.true.)

         clipper = ImGuiListClipper_ImGuiListClipper()
         call ImGuiListClipper_Begin(clipper,n,-1._c_float)
         do while(ImGuiListClipper_Step(clipper))
            call c_f_pointer(clipper,clipper_f)
            do k = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
               call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
               if (cell) then
                  call draw_cp_row(k,f%cpcel(k)%idx,k)
               else
                  call draw_cp_row(k,k,0)
               end if
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call igEndTable()
      end if
    end associate

    ! editing: delete the selected CPs, add one
    call draw_edit_section(w,isys,iview,ttshown)

    ! a nucleus under the mouse: highlight the atom
    if (ihnuc(1) > 0) &
       call sysc(isys)%highlight_atoms(.true.,ihnuc(1:1),ihnuc(2),reshape(ColorHighlightScene,(/4,1/)))

  contains
    !> Row k of the table: symmetry-unique CP i (icp = 0) or cell CP
    !> icp, a copy of symmetry-unique CP i. The properties are those of
    !> the symmetry-unique CP.
    subroutine draw_cp_row(k,i,icp)
      use param, only: bohrtoa
      integer, intent(in) :: k, i, icp

      integer :: j, it
      real*8 :: x(3), xc(3)
      character(len=:), allocatable :: suffix, lbl
      logical :: isbcp, clk, selrow

      associate(c => sys(isys)%c, f => sys(isys)%f(ifield))
        suffix = "_cprow" // string(k)
        it = f%cp(i)%typind
        isbcp = f%isbcp(f%cp(i))

        ! the CP, and the row selectable that highlights it in the view
        ! and, clicked, selects its symmetry-unique CP for deletion (not
        ! the nuclei)
        if (igTableSetColumnIndex(0_c_int)) then
           call iw_text("",alignframe=.true.)
           selrow = .false.
           if (i > c%nneq) selrow = w%cp%sel(i)
           if (iw_highlight_selectable("##cprowsel" // suffix,clicked=clk,selected=selrow)) then
              if (i <= c%nneq) then
                 ! a nucleus: the atom
                 if (icp > 0) then
                    ihnuc = (/icp,atlisttype_ncel_frac/)
                 else
                    ihnuc = (/i,atlisttype_nneq/)
                 end if
              else
                 ihover = (/merge(0,i,icp > 0),icp/)
              end if
           end if
           if (clk .and. i > c%nneq) w%cp%sel(i) = .not.w%cp%sel(i)
           call igSameLine(0._c_float,0._c_float)
           lbl = trim(f%cp(i)%name)
           if (icp > 0) lbl = lbl // " " // string(icp)
           if (i <= c%nneq) then
              if (icp > 0) then
                 call atom_badge(icp,atlisttype_ncel_frac,lbl // "##cprowcp" // suffix)
              else
                 call atom_badge(i,atlisttype_nneq,lbl // "##cprowcp" // suffix)
              end if
           else
              call cp_badge(iview,isys,ifield,i,lbl // "##cprowcp" // suffix)
           end if
        end if

        ! position: fractional for crystals, Cartesian in the input
        ! frame for molecules; the other in the tooltip
        if (icp > 0) then
           x = f%cpcel(icp)%x
        else
           x = f%cp(i)%x
        end if
        xc = (c%x2c(x) + c%molx0) * bohrtoa
        if (igTableSetColumnIndex(1_c_int)) then
           if (ismol) then
              call iw_text(xyz_str(xc))
           else
              call iw_text(xyz_str(x))
              call iw_tooltip("Cartesian: " // xyz_str(xc) // " Å",ttshown)
           end if
        end if

        ! site: multiplicity, Wyckoff letter (crystals), and site symmetry
        if (igTableSetColumnIndex(2_c_int)) then
           lbl = string(f%cp(i)%mult)
           if (allocated(w%cp%wyc) .and. i > c%nneq) lbl = lbl // w%cp%wyc(i-c%nneq)
           call iw_text(lbl // " " // trim(f%cp(i)%pg))
           if (ismol) then
              call iw_tooltip("Multiplicity and site symmetry",ttshown)
           else
              call iw_tooltip("Multiplicity, Wyckoff letter, and site symmetry",ttshown)
           end if
        end if

        ! field, gradient norm, and Laplacian at the CP
        if (igTableSetColumnIndex(3_c_int)) call iw_text(string(f%cp(i)%s%f,'e',decimal=5))
        if (igTableSetColumnIndex(4_c_int)) call iw_text(string(f%cp(i)%s%gfmod,'e',decimal=3))
        if (igTableSetColumnIndex(5_c_int)) call iw_text(string(f%cp(i)%s%del2f,'e',decimal=5))

        ! bond CPs: the ends of the bond path, its length, and the ellipticity
        if (isbcp) then
           do j = 1, 2
              if (igTableSetColumnIndex(int(5+j,c_int))) call path_end_badge(i,icp,j,suffix)
           end do
           if (igTableSetColumnIndex(8_c_int)) then
              if (all(f%cp(i)%ipath > 0)) &
                 call iw_text(string(sum(f%cp(i)%brpathlen) * bohrtoa,'f',decimal=4))
           end if
           if (igTableSetColumnIndex(9_c_int)) then
              if (abs(f%cp(i)%s%hfeval(2)) > 0d0) &
                 call iw_text(string(f%cp(i)%s%hfeval(1)/f%cp(i)%s%hfeval(2)-1d0,'f',decimal=4))
           end if
        end if
      end associate

    end subroutine draw_cp_row

    !> The atom or CP at end j of the bond path of the BCP in row
    !> (symmetry-unique CP i, or cell CP icp).
    subroutine path_end_badge(i,icp,j,suffix)
      integer, intent(in) :: i, icp, j
      character(len=*), intent(in) :: suffix

      integer :: iend, iu
      integer(c_int) :: idx(4)

      associate(c => sys(isys)%c, f => sys(isys)%f(ifield))
        if (icp > 0) then
           ! cell CP: the cell CP at the end, and its lattice vector
           iend = f%cpcel(icp)%ipath(j)
           if (iend >= 1 .and. iend <= c%ncel) then
              ! a nucleus: the cell CP is the cell atom
              idx(1) = iend
              idx(2:4) = f%cpcel(icp)%ilvec(:,j)
              call atom_badge(idx(1),atlisttype_ncel_frac,anchor_label(isys,idx,"?",species=.true.) //&
                 "##cpend" // string(j) // suffix)
              return
           elseif (iend > c%ncel .and. iend <= f%ncpcel) then
              iu = f%cpcel(iend)%idx
              call cp_badge(iview,isys,ifield,iu,trim(f%cp(iu)%name) // " " // string(iend) //&
                 lvec_str(f%cpcel(icp)%ilvec(:,j)) // "##cpend" // string(j) // suffix)
              return
           end if
        else
           ! symmetry-unique CP: the symmetry-unique CP at the end
           iend = f%cp(i)%ipath(j)
           if (iend >= 1 .and. iend <= c%nneq) then
              call atom_badge(iend,atlisttype_nneq,trim(c%at(iend)%name) // "##cpend" // string(j) // suffix)
              return
           elseif (iend > c%nneq .and. iend <= f%ncp) then
              call cp_badge(iview,isys,ifield,iend,trim(f%cp(iend)%name) // "##cpend" // string(j) // suffix)
              return
           end if
        end if

        ! no CP at the end: the unique list tells whether the path
        ! leaves the molecule
        if (f%cp(i)%ipath(j) == -1) then
           call iw_text("(leaves the molecule)",disabled=.true.)
        elseif (f%cp(i)%ipath(j) /= 0) then
           call iw_text("?",disabled=.true.)
        end if
      end associate

    end subroutine path_end_badge

    !> Badge of atom iat of list type itype, in its view color; str is
    !> the label and ImGui ID ("text##id"), as in cp_badge.
    subroutine atom_badge(iat,itype,str)
      integer, intent(in) :: iat, itype
      character(len=*), intent(in) :: str

      real(c_float) :: rgb(3)
      logical :: have, ldum

      have = atom_view_rgb(iview,isys,itype,iat,rgb)
      ldum = iw_atom_button(str,rgb,havergb=have,inert=.true.)

    end subroutine atom_badge

    !> Three coordinates, 4 decimals (no "-0.0000").
    function xyz_str(x) result(str)
      real*8, intent(in) :: x(3)
      character(len=:), allocatable :: str

      real*8 :: xx(3)

      xx = merge(0d0,x,abs(x) < 5d-5)
      str = string(xx(1),'f',decimal=4) // " " // string(xx(2),'f',decimal=4) // " " //&
         string(xx(3),'f',decimal=4)

    end function xyz_str

  end subroutine draw_results_tab

  !> Prepare a run of AUTO or CPREPORT on field ifield of system isys:
  !> connect the system, make the field the reference (plainly:
  !> set_reference would drop the pointprop lists), and use the input
  !> units the options are written in (bohr if the current units have
  !> no length factor, see auto_options). The state changed is saved in
  !> ctx, for auto_leave.
  subroutine auto_enter(isys,ifield,ctx)
    use systemmod, only: sy
    use systems, only: sys
    use global, only: iunit, iunit_bohr, cp_hdegen
    integer, intent(in) :: isys, ifield
    type(auto_context), intent(out) :: ctx

    sy => sys(isys)
    ctx%iref = sys(isys)%iref
    ctx%iunit = iunit
    ctx%hdegen = cp_hdegen
    sys(isys)%iref = ifield
    if (iunit /= 1 .and. iunit /= 2) iunit = iunit_bohr

  end subroutine auto_enter

  !> Restore the state auto_enter changed (AUTO also resets EPSDEGEN).
  subroutine auto_leave(isys,ctx)
    use systems, only: sys
    use global, only: iunit, cp_hdegen
    integer, intent(in) :: isys
    type(auto_context), intent(in) :: ctx

    sys(isys)%iref = ctx%iref
    iunit = ctx%iunit
    cp_hdegen = ctx%hdegen

  end subroutine auto_leave

  !> The Export tab: write the critical points of the field to a file
  !> (CPREPORT): a structure file with the CPs as extra atoms, whose
  !> format is given by the extension, or JSON.
  subroutine draw_export_tab(w,isys,iview,ttshown)
    use systems, only: sys
    use representations, only: field_has_cps
    use utils, only: iw_text, iw_button, iw_tooltip, iw_inputtext, iw_checkbox
    use crystalmod, only: struct_detect_write_format
    use tools_io, only: string, lower, fopen_write, fclose
    use param, only: dirsep, isformat_w_unknown
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(inout) :: ttshown

    integer :: iaux, i0, i1, isformat, lu
    logical :: ldum, isjson
    character(len=:), allocatable :: ext, errexp, line

    if (.not.sys(isys)%goodfield(w%cp%ifield)) then
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       return
    end if
    call iw_text("Field",highlight=.true.)
    call iw_text(string(w%cp%ifield) // ": " // trim(sys(isys)%f(w%cp%ifield)%name),sameline=.true.)
    if (.not.field_has_cps(isys,w%cp%ifield)) then
       call iw_text("This field has no critical points other than the nuclei (search for them in&
          & the Search tab, or load a checkpoint that has them)",disabled=.true.,wrap=.true.)
       return
    end if

    ! file name: editable field plus browse button
    call iw_text("File name",highlight=.true.)
    ldum = iw_inputtext("##cpexportfile",bufsize=1023,texta=w%okfile,width=36)
    call iw_tooltip("File the critical points are written to. The extension gives the format:&
       & a structure file with the critical points as extra atoms (cif, xyz, cri, ...) or JSON (json)",ttshown)
    if (iw_button("Browse...##cpexportbrowse",sameline=.true.)) &
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savecpfile,idparent=w%id,orraise=-1)
    call iw_tooltip("Choose the file with a file browser",ttshown)
    call w%okfile_warn_overwrite()

    ! the extension (json or a structure format)
    i0 = index(w%okfile,dirsep,back=.true.)
    i1 = index(w%okfile(i0+1:),'.',back=.true.)
    ext = ""
    if (i1 > 0) ext = lower(trim(w%okfile(i0+i1+1:)))
    isjson = (ext == "json")

    ! options
    if (.not.isjson) then
       ldum = iw_checkbox("Include the gradient paths (GRAPH)##cpexpgraph",w%cp%expgraph)
       call iw_tooltip("Also write the points of the bond paths, as extra atoms",ttshown)
    end if

    ! write
    errexp = ""
    if (len_trim(w%okfile) == 0) then
       errexp = "Choose a file name"
    elseif (index(trim(w%okfile)," ") > 0) then
       errexp = "The file name cannot contain spaces"
    elseif (len(ext) == 0) then
       errexp = "The file name needs an extension (the format)"
    elseif (.not.isjson) then
       ! a structure format the writer knows (CPREPORT's other
       ! extensions write other files, or a test format)
       call struct_detect_write_format(trim(w%okfile),isformat)
       if (isformat == isformat_w_unknown) errexp = "Unknown file format: ." // ext
    end if
    if (iw_button("Write##cpexportwrite",disabled=(len(errexp) > 0))) then
       ! the writers stop the program if the file cannot be opened:
       ! check first, opening it as they do
       lu = fopen_write(trim(w%okfile),errstop=.false.)
       if (lu < 0) then
          w%errmsg = "Cannot write the file: " // trim(w%okfile)
       else
          call fclose(lu)
          line = trim(w%okfile)
          if (w%cp%expgraph .and. .not.isjson) line = line // " graph"
          call request_job(w,iview,cpjob_export,line)
       end if
    end if
    call iw_tooltip("Write the critical points to the file (the output goes to the output console)",ttshown)
    if (len(errexp) > 0) call iw_text(errexp,danger=.true.,sameline=.true.)

  end subroutine draw_export_tab

  !> Show the seeds of the form of window w in the view of system isys,
  !> computing them again (AUTO, with no search) when the options, the
  !> system, or its geometry changed. The computation waits for the
  !> initialization threads of the system to finish, and for the user
  !> to release the widget being edited (a drag would compute them in
  !> every frame).
  subroutine show_seed_preview(w,isys)
    use autocp, only: autocritic
    use systems, only: sys, sysc, sys_ready, ok_system, are_threads_running
    use representations, only: rep_shape, shapekind_sphere
    use global, only: iunit
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    integer :: i, n
    logical :: ok, stale, found, newseeds
    character(len=:), allocatable :: line
    type(rep_shape), allocatable :: shp(:)
    type(auto_context) :: ctx

    ! the options (empty if the form cannot make valid ones)
    line = ""
    if (size(w%cp%seed) > 0) then
       if (len(form_error(w)) == 0) line = auto_options(w,isys,iunit)
    end if

    ! compute the seeds again if something changed (absolute frame)
    stale = .not.allocated(w%cp%seedline)
    if (.not.stale) stale = (w%cp%seedsys /= isys) .or. (w%cp%seedline /= line) .or.&
       (w%cp%seedtime /= sysc(isys)%timelastchange_geometry)
    newseeds = .false.
    if (stale .and. ok_system(isys,sys_ready) .and. .not.are_threads_running() .and.&
       .not.igIsAnyItemActive()) then
       if (allocated(w%cp%seedx)) deallocate(w%cp%seedx)
       if (len(line) > 0) then
          call auto_enter(isys,sys(isys)%iref,ctx)
          call autocritic(line,ok,seeds=w%cp%seedx)
          call auto_leave(isys,ctx)
          if (.not.ok .and. allocated(w%cp%seedx)) deallocate(w%cp%seedx)
          if (allocated(w%cp%seedx)) w%cp%seedx = w%cp%seedx + spread(sys(isys)%c%molx0,2,size(w%cp%seedx,2))
       end if
       w%cp%seedline = line
       w%cp%seedsys = isys
       w%cp%seedtime = sysc(isys)%timelastchange_geometry
       newseeds = .true.
    end if

    ! draw the last seeds computed for this system: keep the shapes in
    ! the view, and build them again only if they changed or are gone
    if (.not.allocated(w%cp%seedx) .or. w%cp%seedsys /= isys) return
    if (sysc(isys)%sc%isinit == 0) return
    n = min(size(w%cp%seedx,2),maxseedshow)
    if (n == 0) return
    call sysc(isys)%sc%show_transient_shapes(w%id,1,found=found)
    if (found .and. .not.newseeds) return
    allocate(shp(n))
    do i = 1, n
       shp(i) = rep_shape(kind=shapekind_sphere,x1=w%cp%seedx(:,i),rad=seed_rad,rgb=seed_rgb)
    end do
    call sysc(isys)%sc%show_transient_shapes(w%id,1,shp)

  end subroutine show_seed_preview

  !> The editing of the CP list, under the results table: delete the
  !> selected CPs, or add one by a search from a point (the center of
  !> the selected atoms, or a bond or point picked in the view).
  subroutine draw_edit_section(w,isys,iview,ttshown)
    use systems, only: sys, sysc
    use utils, only: iw_text, iw_button, iw_tooltip
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(inout) :: ttshown

    integer :: nsel, nat, k
    integer, allocatable :: iat(:)
    real*8 :: x1(3), xd(3), xs(3)
    character(len=:), allocatable :: errform

    ! delete
    nsel = count(w%cp%sel)
    if (iw_button("Delete selected##cpdelete",disabled=(nsel == 0))) &
       call request_job(w,iview,cpjob_delete)
    call iw_tooltip("Delete the selected critical points (click the rows to select them; the nuclei&
       & cannot be deleted) with all their copies in the cell, and trace the bond paths again",ttshown)
    if (iw_button("Clear selection##cpclearsel",sameline=.true.,disabled=(nsel == 0))) w%cp%sel = .false.
    call iw_tooltip("Unselect all critical points",ttshown)
    if (nsel > 0) call iw_text(string(nsel) // " selected",sameline=.true.)

    ! add: a search from one point, with the advanced options of the
    ! Search tab
    call iw_text("Add a critical point",highlight=.true.)
    call iw_tooltip("Search for a critical point from one point, with the advanced options of&
       & the Search tab. The new critical point is added to the list",ttshown)
    errform = form_error(w)
    if (iw_button("From selection##cpaddsel",disabled=(len(errform) > 0 .or. w%cp%picking))) then
       call sysc(isys)%highlighted_atom_list(nat,iat)
       if (nat == 0) then
          w%errmsg = "Select some atoms first: the search starts at their center"
       else
          ! the center of the selected atoms, each taken at its image
          ! nearest the first one (crystals)
          associate(c => sys(isys)%c)
            x1 = c%atcel(iat(1))%x
            xs = 0d0
            do k = 2, nat
               xd = c%atcel(iat(k))%x - x1
               if (.not.c%ismolecule) xd = xd - nint(xd)
               xs = xs + xd
            end do
            call request_add(w,isys,iview,c%x2c(x1 + xs / nat))
          end associate
       end if
    end if
    call iw_tooltip("Search for a critical point starting at the center of the selected atoms",ttshown)
    if (iw_button("Pick in view##cpaddpick",sameline=.true.,disabled=(len(errform) > 0 .or. w%cp%picking))) then
       w%cp%picking = .true.
       w%cp%pickview = iview
       call w%cp%pick%arm()
       call win(iview)%viewmode_set_forced(vm_pick_bond,"Pick a bond or a point to search for a&
          & critical point from",w%id,acceptempty=.true.)
    end if
    call iw_tooltip("Search for a critical point starting at a point picked in the view: the middle&
       & of a bond if one is clicked, otherwise the clicked point (on the plane through the center of&
       & the scene)",ttshown)
    if (w%cp%picking) call iw_text("Click a bond or a point in the view",disabled=.true.,sameline=.true.)
    if (len(errform) > 0) call iw_text(errform // " (Search tab, advanced options)",danger=.true.,wrap=.true.)

  end subroutine draw_edit_section

  !> Handle the pending pick of the point to add a CP from, commanded
  !> to view iview (system isys): a picked bond gives its midpoint,
  !> anything else the clicked point. Then the add job is requested.
  subroutine poll_add_pick(w,isys,iview)
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview

    integer :: istat
    real*8 :: xc(3)

    call view_pick_result(iview,w%id,isys,w%cp%pick,istat,xc)
    if (istat == ipick_point) call request_add(w,isys,iview,xc)
    if (istat /= ipick_pending) w%cp%picking = .false.

  end subroutine poll_add_pick

  !> Request the job that adds a CP to the list by a search from xc
  !> (cell-frame Cartesian bohr), for the view iview of system isys.
  subroutine request_add(w,isys,iview,xc)
    use systems, only: sys
    use global, only: iunit
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    real*8, intent(in) :: xc(3)

    real*8 :: x(3)

    ! in the coordinates of the form
    associate(c => sys(isys)%c)
      if (c%ismolecule) then
         x = (xc + c%molx0) * bohrtoa
      else
         x = c%c2x(xc)
      end if
    end associate
    call request_job(w,iview,cpjob_add,auto_options(w,isys,iunit,point=x))

  end subroutine request_add

  !> Request a blocking job of kind kind (cpjob_*) for view iview, with
  !> the AUTO options line (search, add). The main loop runs it after
  !> the frame with its overlay; a second request in the same frame is
  !> ignored.
  subroutine request_job(w,iview,kind,line)
    use gui_main, only: pending_block_window
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview, kind
    character(len=*), intent(in), optional :: line

    if (pending_block_window /= 0) return
    w%cp%pending_kind = kind
    if (present(line)) w%cp%pending_line = line
    w%cp%pending_view = iview
    w%errmsg = ""
    pending_block_window = w%id

  end subroutine request_job

  !> Recompute the caches of the results table if the CP list of the
  !> field (or the field) changed: the Wyckoff letters of the
  !> symmetry-unique CPs and the summary. The selection is cleared.
  subroutine update_table_caches(w,isys)
    use systems, only: sys, sysc
    use representations, only: cp_wyckoff
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    if (w%cp%tfield == w%cp%ifield .and. w%cp%ttime == sysc(isys)%timelastchange_cplist) return
    w%cp%tfield = w%cp%ifield
    w%cp%ttime = sysc(isys)%timelastchange_cplist
    call cp_wyckoff(isys,w%cp%ifield,w%cp%wyc)
    w%cp%summary = cp_summary(isys,w%cp%ifield)
    if (allocated(w%cp%sel)) deallocate(w%cp%sel)
    allocate(w%cp%sel(sys(isys)%f(w%cp%ifield)%ncp))
    w%cp%sel = .false.

  end subroutine update_table_caches

  !> Highlight in the critical points object of view iview (for the
  !> window's field) the CP ihover (as rep_cps%ihover), and stop
  !> highlighting the previous one. Only changes rebuild the lists.
  subroutine set_cp_hover(w,iview,ihover)
    use representations, only: reptype_cps
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview
    integer, intent(in) :: ihover(2)

    if (all(ihover == w%cp%ihover) .and. iview == w%cp%hoverview) return
    call apply(w%cp%hoverview,(/0,0/))
    call apply(iview,ihover)
    w%cp%ihover = ihover
    w%cp%hoverview = iview

  contains
    subroutine apply(iv,ih)
      integer, intent(in) :: iv, ih(2)

      integer :: i
      logical :: ch

      if (iv < 1 .or. iv > nwin) return
      if (.not.win(iv)%isopen .or. .not.associated(win(iv)%sc)) return
      ch = .false.
      associate(sc => win(iv)%sc)
        do i = 1, sc%nrep
           if (.not.sc%rep(i)%isinit .or. sc%rep(i)%type /= reptype_cps) cycle
           ! set it in the object of the field; clear it in all
           if (any(ih /= 0) .and. sc%rep(i)%cps%ifield /= w%cp%ifield) cycle
           if (all(sc%rep(i)%cps%ihover == ih)) cycle
           sc%rep(i)%cps%ihover = ih
           ch = .true.
        end do
        if (ch) sc%forcebuildlists = .true.
      end associate

    end subroutine apply
  end subroutine set_cp_hover

  !> The Search tab: field, seeds, advanced options, and the Run button.
  subroutine draw_search_tab(w,isys,iview,ttshown)
    use systems, only: sys
    use global, only: iunit
    use gui_main, only: g
    use tools_io, only: string
    use utils, only: iw_text, iw_button, iw_tooltip, iw_field_combo, iw_calcwidth,&
       iw_checkbox, iw_helpermark
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(inout) :: ttshown

    integer :: i, idel
    integer(c_int) :: ifield
    logical :: ldum
    real(c_float) :: xcol
    character(len=:), allocatable :: xunit, errrun
    character(kind=c_char,len=:), allocatable, target :: str1
    type(ImVec2) :: sz

    ! maximum height of the seed list, in rows
    integer, parameter :: maxrow_seeds = 12

    if (sys(isys)%c%ismolecule) then
       xunit = " (Å)"
    else
       xunit = " (fractional)"
    end if

    ! field
    call iw_text("Field",highlight=.true.,alignframe=.true.)
    call igSameLine(0._c_float,-1._c_float)
    ifield = w%cp%ifield
    if (iw_field_combo("##cpfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
       nonestr="<field not available>")) w%cp%ifield = ifield
    call iw_tooltip("Field whose critical points are searched",ttshown)
    if (.not.sys(isys)%goodfield(w%cp%ifield)) then
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)
       return
    end if

    ! the seeds, in a scrolling box of at most maxrow_seeds rows
    call iw_text("Seeds",highlight=.true.)
    call iw_helpermark("The search for critical points starts from a list of points, the&
       & seeds. Each entry in this list adds a set of seeds. A Newton-Raphson search runs from every seed, and&
       & the point it converges to is added to the list of critical points.",sameline=.true.)

    ! the height of the box is that of its content in the last frame
    ! (the bundled imgui cannot fit a child to its content), from one
    ! row to maxrow_seeds rows
    sz%x = 0._c_float
    sz%y = max(w%cp%seedh,igGetFrameHeightWithSpacing() + 2 * g%Style%WindowPadding%y)
    sz%y = min(sz%y,maxrow_seeds * igGetFrameHeightWithSpacing() + 2 * g%Style%WindowPadding%y)
    str1 = "##cpseedlist" // c_null_char
    idel = 0
    if (igBeginChild_Str(c_loc(str1),sz,.true._c_bool,ImGuiWindowFlags_None)) then
       ! the parameter widgets start at one column, after the longest label
       xcol = igGetCursorPosX() + g%Style%IndentSpacing + seed_label_width(xunit) + 2 * g%Style%ItemSpacing%x
       do i = 1, size(w%cp%seed)
          if (i > 1) call igSeparator()
          if (draw_seed(w,isys,i,xunit,xcol,ttshown)) idel = i
       end do
       w%cp%seedh = igGetCursorPosY() - g%Style%ItemSpacing%y + g%Style%WindowPadding%y
    end if
    call igEndChild()
    if (idel > 0) w%cp%seed = [w%cp%seed(:idel-1), w%cp%seed(idel+1:)]
    if (iw_button("Add seed##cpseedadd")) call add_seed(w,isys,cpseed_point)
    call iw_tooltip("Add a seed (a point at the center of the system; change its kind above)",ttshown)
    if (iw_button("Reset seeds##cpseedreset",sameline=.true.)) call reset_seeds(w,isys)
    call iw_tooltip("Go back to the seeding AUTO uses by default for this system",ttshown)
    ldum = iw_checkbox("Show seeds##cpshowseeds",w%cp%showseeds)
    call iw_tooltip("Show the starting points of the searches in the view (grey spheres), after&
       & the pruning and the clipping",ttshown)
    if (w%cp%showseeds .and. allocated(w%cp%seedx) .and. w%cp%seedsys == isys) then
       if (size(w%cp%seedx,2) > maxseedshow) then
          call iw_text("(" // string(size(w%cp%seedx,2)) // " seeds, showing " // string(maxseedshow) // ")",&
             sameline=.true.)
       else
          call iw_text("(" // string(size(w%cp%seedx,2)) // " seeds)",sameline=.true.)
       end if
    end if

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
    if (size(w%cp%seed) == 0) then
       errrun = "Add a seed"
    else
       errrun = form_error(w)
    end if
    if (iw_button("Run##cprun",danger=.true.,disabled=(len(errrun) > 0))) &
       call request_job(w,iview,cpjob_search,auto_options(w,isys,iunit))
    call iw_tooltip("Search for the critical points (the calculation blocks the interface&
       & until it finishes; the output goes to the output console)",ttshown)
    if (len(errrun) > 0) call iw_text(errrun,danger=.true.,sameline=.true.)

    ! the result of the last run
    if (allocated(w%cp%summary)) call iw_text(w%cp%summary,wrap=.true.)

  end subroutine draw_search_tab

  !> Width of the longest label of the seed parameters (xunit = the
  !> units of the positions).
  function seed_label_width(xunit) result(wid)
    use utils, only: iw_textwidth
    character(len=*), intent(in) :: xunit
    real(c_float) :: wid

    wid = max(iw_textwidth("Subdivision level"),iw_textwidth("Center" // xunit),&
       iw_textwidth("Radius (Å)"),iw_textwidth("Radial points"),iw_textwidth("Polar points"),&
       iw_textwidth("Azimuthal points"),iw_textwidth("Maximum distance (Å)"),&
       iw_textwidth("Points per pair"),iw_textwidth("Start" // xunit),iw_textwidth("End" // xunit),&
       iw_textwidth("Points"),iw_textwidth("Position" // xunit))

  end function seed_label_width

  !> Draw the widgets of seed i of window w: its kind, a button to
  !> remove it (returns true if pressed), and the parameters that kind
  !> uses, each with its label on the left and the widget at column
  !> xcol. Choosing a kind gives its parameters usable values.
  function draw_seed(w,isys,i,xunit,xcol,ttshown) result(del)
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_inputint, iw_dragfloat_real8,&
       iw_close_button
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, i
    character(len=*), intent(in) :: xunit
    real(c_float), intent(in) :: xcol
    logical, intent(inout) :: ttshown
    logical :: del

    character(len=:), allocatable :: suf, tt
    logical :: ldum, ch

    suf = "##cpseed" // string(i)
    call iw_text(string(i) // ".",alignframe=.true.)
    call iw_combo_simple("##cpseedkind" // suf,seedkind_names,w%cp%seed(i)%typ,sameline=.true.,&
       changed=ch,tooltips=seedkind_tooltips,ttshown=ttshown)
    if (ch) call kind_defaults(w,isys,i)
    call igSameLine(0._c_float,-1._c_float)
    del = iw_close_button("##cpseeddel" // suf)
    call iw_tooltip("Remove this seed",ttshown)

    call igIndent(0._c_float)
    associate(s => w%cp%seed(i))
      select case(s%typ)
      case (cpseed_ws,cpseed_oh)
         tt = "Number of recursive subdivisions of the region (0 to 7; more seeds at each level)"
         call label("Subdivision level",tt)
         ldum = iw_inputint(suf // "dep",s%depth,width=8)
         s%depth = max(0,min(s%depth,maxdepth))
         call iw_tooltip(tt,ttshown)
         call label("Center" // xunit)
         call x0_widget()
         if (s%typ == cpseed_ws) then
            tt = "Scale the Wigner-Seitz cell to this radius around the center (0 = the whole cell)"
            call label("Radius (Å)",tt)
            ldum = iw_dragfloat_real8(suf // "rad",x1=s%rad,speed=0.01d0,min=0d0,max=100d0,decimal=3)
            call iw_tooltip(tt,ttshown)
            if (s%rad <= 0d0) call iw_text("(the whole cell)",disabled=.true.,sameline=.true.)
         else
            call label("Radius (Å)")
            ldum = iw_dragfloat_real8(suf // "rad",x1=s%rad,speed=0.01d0,min=0.01d0,max=100d0,decimal=3)
            call label("Radial points")
            ldum = iw_inputint(suf // "nr",s%nr,width=8)
            s%nr = max(s%nr,1)
         end if
      case (cpseed_sphere)
         call label("Center" // xunit)
         call x0_widget()
         call label("Radius (Å)")
         ldum = iw_dragfloat_real8(suf // "rad",x1=s%rad,speed=0.01d0,min=0.01d0,max=100d0,decimal=3)
         call label("Polar points")
         ldum = iw_inputint(suf // "nth",s%ntheta,width=8)
         call label("Azimuthal points")
         ldum = iw_inputint(suf // "nph",s%nphi,width=8)
         call label("Radial points")
         ldum = iw_inputint(suf // "nr",s%nr,width=8)
         s%ntheta = max(s%ntheta,1)
         s%nphi = max(s%nphi,1)
         s%nr = max(s%nr,1)
      case (cpseed_pair,cpseed_triplet)
         tt = "Only atoms closer than this distance are combined"
         call label("Maximum distance (Å)",tt)
         ldum = iw_dragfloat_real8(suf // "dist",x1=s%dist,speed=0.01d0,min=0d0,max=100d0,decimal=3)
         call iw_tooltip(tt,ttshown)
         if (s%typ == cpseed_pair) then
            tt = "Number of seeds on the segment between the two atoms"
            call label("Points per pair",tt)
            ldum = iw_inputint(suf // "npts",s%npts,width=8)
            s%npts = max(s%npts,1)
            call iw_tooltip(tt,ttshown)
         end if
      case (cpseed_line)
         call label("Start" // xunit)
         call x0_widget()
         call label("End" // xunit)
         ldum = iw_dragfloat_real8(suf // "x1",x3=s%x1,speed=0.001d0,decimal=4)
         call label("Points")
         ldum = iw_inputint(suf // "npts",s%npts,width=8)
         s%npts = max(s%npts,2)
      case (cpseed_point)
         call label("Position" // xunit)
         call x0_widget()
      case (cpseed_mesh)
         call iw_text("The points of the molecular integration mesh",disabled=.true.)
      end select
    end associate
    call igUnindent(0._c_float)

  contains
    !> A label in the label column, with tooltip tt; the widget follows
    !> at xcol.
    subroutine label(str,tt)
      character(len=*), intent(in) :: str
      character(len=*), intent(in), optional :: tt
      call iw_text(str,alignframe=.true.)
      ! (no delay: text has no ID, so the delayed tooltip would not show)
      if (present(tt)) call iw_tooltip(tt)
      call igSameLine(xcol,-1._c_float)
    end subroutine label
    subroutine x0_widget()
      ldum = iw_dragfloat_real8(suf // "x0",x3=w%cp%seed(i)%x0,speed=0.001d0,decimal=4)
    end subroutine x0_widget
  end function draw_seed

  !> The advanced options of AUTO, in a collapsed tree node; each one
  !> is used only if its checkbox is set (otherwise, AUTO's default).
  !> The labels are in one column (after the checkbox, if any), and the
  !> widgets start at another.
  subroutine draw_advanced(w,isys,xunit,ttshown)
    use gui_main, only: g
    use utils, only: iw_tooltip, iw_checkbox, iw_inputtext, iw_dragfloat_real8,&
       iw_combo_simple, iw_text, iw_textwidth
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    character(len=*), intent(in) :: xunit
    logical, intent(inout) :: ttshown

    character(len=:,kind=c_char), allocatable, target :: strad
    character(len=:), allocatable :: tt
    logical :: ldum, ch
    real(c_float) :: xlab, xcol

    strad = "Advanced options##cpadvanced" // c_null_char
    if (igTreeNodeEx_Str(c_loc(strad),ImGuiTreeNodeFlags_None)) then
       ! the label column (after a checkbox) and the widget column
       xlab = igGetCursorPosX() + igGetFrameHeight() + g%Style%ItemInnerSpacing%x
       xcol = xlab + 2 * g%Style%ItemSpacing%x + max(iw_textwidth("Gradient norm (GRADEPS)"),&
          iw_textwidth("CP distance (CPEPS, Å)"),iw_textwidth("Nucleus distance (NUCEPS, Å)"),&
          iw_textwidth("Hydrogen distance (NUCEPSH, Å)"),iw_textwidth("Degenerate eigenvalue (EPSDEGEN)"),&
          iw_textwidth("Discard (DISCARD)"),iw_textwidth("Clip (CLIP)"),iw_textwidth("Box corner 1" // xunit),&
          iw_textwidth("Center" // xunit),iw_textwidth("Radius (Å)"))

       tt = "Maximum gradient norm of a critical point (default 1e-12)"
       call check_label("Gradient norm (GRADEPS)##cpusegradeps",w%cp%use_gradeps,tt)
       ldum = iw_inputtext("##cpgradeps",bufsize=31,textf=w%cp%gradeps,width=8)
       call iw_tooltip(tt,ttshown)

       tt = "Two critical points closer than this are the same (default 0.0053 Å)"
       call check_label("CP distance (CPEPS, Å)##cpusecpeps",w%cp%use_cpeps,tt)
       ldum = iw_dragfloat_real8("##cpcpeps",x1=w%cp%cpeps,speed=0.001d0,min=0d0,max=10d0,decimal=4)
       call iw_tooltip(tt,ttshown)

       tt = "Discard critical points closer than this to a nucleus (default 0.053 Å; larger for grids)"
       call check_label("Nucleus distance (NUCEPS, Å)##cpusenuceps",w%cp%use_nuceps,tt)
       ldum = iw_dragfloat_real8("##cpnuceps",x1=w%cp%nuceps,speed=0.001d0,min=0d0,max=10d0,decimal=4)
       call iw_tooltip(tt,ttshown)

       tt = "Same, for hydrogen nuclei (default 0.106 Å; larger for grids)"
       call check_label("Hydrogen distance (NUCEPSH, Å)##cpusenucepsh",w%cp%use_nucepsh,tt)
       ldum = iw_dragfloat_real8("##cpnucepsh",x1=w%cp%nucepsh,speed=0.001d0,min=0d0,max=10d0,decimal=4)
       call iw_tooltip(tt,ttshown)

       tt = "Hessian eigenvalues smaller than this make a degenerate critical point (default 1e-8)"
       call check_label("Degenerate eigenvalue (EPSDEGEN)##cpuseepsdegen",w%cp%use_epsdegen,tt)
       ldum = iw_inputtext("##cpepsdegen",bufsize=31,textf=w%cp%epsdegen,width=8)
       call iw_tooltip(tt,ttshown)

       tt = "Discard the critical points where this expression is not zero, e.g. $0 < 1e-3 (empty = none)"
       call label("Discard (DISCARD)",tt)
       ldum = iw_inputtext("##cpdiscard",bufsize=1023,textf=w%cp%discard,width=30)
       call iw_tooltip(tt,ttshown)

       tt = "Only use the seeds inside this region"
       call label("Clip (CLIP)",tt)
       call iw_combo_simple("##cpclip","None" // c_null_char // "Box" // c_null_char //&
          "Sphere" // c_null_char,w%cp%iclip,changed=ch)
       call iw_tooltip(tt,ttshown)
       if (ch) call clip_defaults(w,isys)
       if (w%cp%iclip == 1) then
          call label("Box corner 1" // xunit)
          ldum = iw_dragfloat_real8("##cpclipx0",x3=w%cp%clipx0,speed=0.001d0,decimal=4)
          call label("Box corner 2" // xunit)
          ldum = iw_dragfloat_real8("##cpclipx1",x3=w%cp%clipx1,speed=0.001d0,decimal=4)
       elseif (w%cp%iclip == 2) then
          call label("Center" // xunit)
          ldum = iw_dragfloat_real8("##cpclipx0s",x3=w%cp%clipx0,speed=0.001d0,decimal=4)
          call label("Radius (Å)")
          ldum = iw_dragfloat_real8("##cpcliprad",x1=w%cp%cliprad,speed=0.01d0,min=0d0,max=100d0,decimal=3)
       end if
       call igTreePop()
    end if

  contains
    !> A checkbox carrying the label (str) of an option that is used
    !> only if checked, with tooltip tt; the widget follows at xcol.
    subroutine check_label(str,flag,tt)
      character(len=*), intent(in) :: str
      logical, intent(inout) :: flag
      character(len=*), intent(in) :: tt
      logical :: ld
      ld = iw_checkbox(str,flag)
      call iw_tooltip(tt,ttshown)
      call igSameLine(xcol,-1._c_float)
    end subroutine check_label
    !> A label in the label column, without a checkbox, with tooltip tt
    !> if given; the widget follows at xcol.
    subroutine label(str,tt)
      character(len=*), intent(in) :: str
      character(len=*), intent(in), optional :: tt
      call igSetCursorPosX(xlab)
      call iw_text(str,alignframe=.true.)
      ! (no delay: text has no ID, so the delayed tooltip would not show)
      if (present(tt)) call iw_tooltip(tt)
      call igSameLine(xcol,-1._c_float)
    end subroutine label
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

  !> Why the advanced options of window w would not make valid AUTO
  !> options, or an empty string if they do.
  function form_error(w) result(errmsg)
    type(window), intent(in) :: w
    character(len=:), allocatable :: errmsg

    errmsg = ""
    if (w%cp%use_gradeps .and. .not.isnumber(w%cp%gradeps)) then
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
  !> in bohr, as AUTO reads them. If point is present (a position in
  !> the coordinates of the form), the only seed is a POINT there, and
  !> there is no CLIP (to add a CP from that point).
  function auto_options(w,isys,iunit,point) result(line)
    use systems, only: sys
    use global, only: dunit0
    use tools_io, only: string
    type(window), intent(in) :: w
    integer, intent(in) :: isys, iunit
    real*8, intent(in), optional :: point(3)
    character(len=:), allocatable :: line

    logical :: ismol
    real*8 :: fl

    ismol = sys(isys)%c%ismolecule
    fl = 1d0 / bohrtoa ! Å to the input units
    if (iunit == 1 .or. iunit == 2) fl = fl * dunit0(iunit)

    ! seeds: the point, or those of the form
    line = ""
    if (present(point)) then
       line = " seed point" // xstr(" x0",point)
    else
       call form_seeds()
    end if

    ! advanced options (no CLIP for a point)
    if (w%cp%use_gradeps) line = line // " gradeps " // trim(adjustl(w%cp%gradeps))
    if (w%cp%use_cpeps) line = line // " cpeps " // rstr(w%cp%cpeps / bohrtoa)
    if (w%cp%use_nuceps) line = line // " nuceps " // rstr(w%cp%nuceps / bohrtoa)
    if (w%cp%use_nucepsh) line = line // " nucepsh " // rstr(w%cp%nucepsh / bohrtoa)
    if (w%cp%use_epsdegen) line = line // " epsdegen " // trim(adjustl(w%cp%epsdegen))
    if (len_trim(w%cp%discard) > 0) line = line // ' discard "' // trim(w%cp%discard) // '"'
    if (.not.present(point)) then
       if (w%cp%iclip == 1) then
          line = line // " clip cube" // xstr("",w%cp%clipx0) // xstr("",w%cp%clipx1)
       elseif (w%cp%iclip == 2) then
          line = line // " clip sphere" // xstr("",w%cp%clipx0) // " " // rstr(w%cp%cliprad * fl)
       end if
    end if
    line = line(2:)

  contains
    !> Append the seeds of the form to line.
    subroutine form_seeds()
      integer :: i

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

    end subroutine form_seeds

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
