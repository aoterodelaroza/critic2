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
  use meshmod, only: mesh_level_kw
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
     "The points of the molecular integration mesh of the symmetry-unique atoms, each atom's grid&
     & pruned to its Voronoi region and without its core." //&
     c_null_char

  ! the format combo of the Export tab: auto-detect (0), JSON, a VMD
  ! script, then the structure formats in fmtperm order (built on first use)
  integer, parameter :: iexp_json = 1
  integer, parameter :: iexp_vmd = 2
  integer, parameter :: iexp_fmt1 = 3
  character(len=17), parameter :: expformat_names(2) = (/"JSON file (.json)","VMD script (.vmd)"/)
  character(len=4), parameter :: expformat_exts(2) = (/"json","vmd "/)
  character(kind=c_char,len=:), allocatable :: expformat_combostr

  ! levels of the mesh of a MESH seed (the order of mlevel and
  ! mesh_level_kw), and their descriptions
  character(len=*), parameter :: meshlevel_names = "Small" // c_null_char // "Normal" // c_null_char //&
     "Good" // c_null_char // "Very good" // c_null_char // "Amazing" // c_null_char
  character(len=*), parameter :: meshlevel_tooltips = &
     "The coarsest atomic grids: the fewest seeds (the default)" // c_null_char //&
     "Finer atomic grids: about twice the seeds of Small" // c_null_char //&
     "The atomic grids of integrations: about four times the seeds of Small" // c_null_char //&
     "Very fine atomic grids: many seeds" // c_null_char //&
     "The finest atomic grids: very many seeds" // c_null_char

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
          call cancel_pick(w)
          ! the default export file: <root>.cps.cif (<root>.cif could be the source file)
          w%okfile = okfile_default(isys,"structure","cps.cif")
          w%cp%tfield = -1
       end if
    end if

    ! a pending pick in the view it was armed on (the point to add a CP
    ! from, or a position of the form); dropped, releasing that view, if
    ! the window closes, the system is not ready, or the window moved to
    ! another view
    if (w%cp%picking > 0) then
       if (goodsys .and. .not.doquit .and. iview == w%cp%pickview) then
          call poll_pick(w,isys,iview)
       else
          call cancel_pick(w)
       end if
    end if

    ! the system
    ihover = 0
    if (goodsys) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! the field, shared by the three tabs (they draw nothing more if
       ! it is not available)
       call draw_field_combo(w,isys,ttshown)

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
    if (goodsys .and. .not.doquit .and. w%cp%showseeds) then
       call update_seeds(w,isys)
       call draw_seed_preview(w,isys,iview)
    end if

    ! the clip region of the form, shown in the view
    if (goodsys .and. .not.doquit .and. w%cp%iclip > 0) call draw_clip_preview(w,isys,iview)

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
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_wait_overlay
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    character(len=:), allocatable :: title, info
    integer :: isys

    ! what is being done
    select case (w%cp%pending_kind)
    case (cpjob_add)
       title = "Searching for a critical point..."
    case (cpjob_delete)
       title = "Deleting critical points and tracing the bond paths..."
    case (cpjob_export)
       title = "Writing the critical points..."
    case (cpjob_estimate)
       title = "Estimating the time of the search..."
    case default
       title = "Searching for critical points..."
    end select

    ! the system, field, and input
    info = ""
    isys = w%cp%isys
    if (ok_system(isys,sys_init)) then
       info = "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (sys(isys)%goodfield(w%cp%ifield)) &
          info = info // newline // "Field:  " // string(w%cp%ifield) // ": " //&
          trim(sys(isys)%f(w%cp%ifield)%name)
    end if
    if (allocated(w%cp%pending_line)) then
       select case (w%cp%pending_kind)
       case (cpjob_search,cpjob_add)
          info = info // newline // "Input:  AUTO " // w%cp%pending_line
       case (cpjob_export)
          info = info // newline // "Input:  CPREPORT " // w%cp%pending_line
       end select
    end if

    call iw_wait_overlay(title,info,cp_cancellable(w%cp%pending_kind))

  end subroutine block_cp

  !> Whether a job of this kind (cpjob_*) can be cancelled with Esc:
  !> the searches and the deletion (whose bond paths are traced
  !> again).
  logical function cp_cancellable(kind)
    integer, intent(in) :: kind

    cp_cancellable = (kind == cpjob_search .or. kind == cpjob_delete)

  end function cp_cancellable

  !> Run the blocking job of the critical points window, called by the
  !> main loop after the frame with the overlay, on the chosen field: a
  !> search (AUTO with the options of the Search tab), the addition of
  !> a CP (AUTO from one point, appending), or the deletion of the CPs
  !> in pending_del (and the bond graph traced again). Then the
  !> checkpoint and the critical points and gradient paths objects in
  !> the view. If AUTO rejects the options, the field
  !> is left as it was.
  module subroutine run_cp_pending(w)
    use autocp, only: autocritic, autocritic_graph, cpreport
    use gui_main, only: begin_cancellable, end_cancellable
    use systems, only: sys, sysc, sys_init, ok_system, lastchange_cplist
    use tools_io, only: uout, string
    use fieldmod, only: cplist_backup
    class(window), intent(inout), target :: w

    integer :: isys, ifield, iview, ncp0, kind, lu, ios
    type(auto_context) :: ctx
    type(cplist_backup) :: cpback
    logical :: ok, changes, cancelled, lex
    character(len=:), allocatable :: cpfile, errmsg

    isys = w%cp%isys
    ifield = w%cp%ifield
    iview = w%cp%pending_view
    kind = w%cp%pending_kind
    ok = ok_system(isys,sys_init)
    if (ok) ok = sys(isys)%goodfield(ifield)
    if (ok) then
       if (kind == cpjob_delete) then
          ok = allocated(w%cp%pending_del)
          if (ok) ok = (size(w%cp%pending_del) == sys(isys)%f(ifield)%ncp)
       elseif (kind == cpjob_estimate) then
          ok = allocated(w%cp%seedx) .and. (w%cp%seedsys == isys) .and. (w%cp%seedfield == ifield)
       else
          ok = allocated(w%cp%pending_line)
       end if
    end if
    if (.not.ok) then
       w%errmsg = "The system or the field is no longer available"
       return
    end if

    ! run on the chosen field (the main loop stops the initialization
    ! threads and reads the output)
    call auto_enter(isys,ifield,ctx)
    cancelled = .false.
    ncp0 = sys(isys)%f(ifield)%ncp
    if (cp_cancellable(kind)) call begin_cancellable()
    if (kind == cpjob_delete) then
       ! a cancel restores the CP list (AUTO restores it itself)
       write (uout,'("* Deleting ",A," critical points (GUI) and tracing the bond paths again")') &
          string(count(w%cp%pending_del))
       call sys(isys)%f(ifield)%backup_cplist(cpback)
       call sys(isys)%f(ifield)%delete_cps(w%cp%pending_del)
       call autocritic_graph()
       ok = .true.
    elseif (kind == cpjob_export) then
       call cpreport(w%cp%pending_line)
       write (uout,'("* Critical points written to: ",A/)') trim(w%okfile)
       ok = .true.
    elseif (kind == cpjob_estimate) then
       call estimate_search(w,isys,ifield)
       ok = .true.
    else
       call autocritic(w%cp%pending_line,ok,clear=(kind == cpjob_search .and. w%cp%discard_existing))
    end if
    if (cp_cancellable(kind)) cancelled = end_cancellable()
    if (cancelled .and. kind == cpjob_delete) then
       call sys(isys)%f(ifield)%restore_cplist(cpback)
       write (uout,'("+ The deletion was cancelled: the list of critical points is unchanged"/)')
    elseif (kind == cpjob_delete) then
       write (uout,*)
    end if
    if (cancelled) ok = .false.
    call auto_leave(isys,ctx)

    ! the checkpoint, unless disabled (not AUTO's CHK, which reads an
    ! existing checkpoint instead of searching)
    changes = (kind /= cpjob_export .and. kind /= cpjob_estimate)
    if (ok .and. .not.w%cp%nochk .and. changes) then
       cpfile = sys(isys)%f(ifield)%chk_cps_file()
       if (len(cpfile) == 0) then
          continue
       elseif (sys(isys)%f(ifield)%ncp <= sys(isys)%c%nneq) then
          ! only the nuclei are left: remove the checkpoint instead
          inquire(file=cpfile,exist=lex)
          if (lex) then
             open(newunit=lu,file=cpfile,status="old",iostat=ios)
             if (ios == 0) close(lu,status="delete",iostat=ios)
             if (ios == 0) then
                write (uout,'("* Checkpoint file removed: ",A/)') cpfile
             else
                write (uout,'("!! Warning !! Could not remove the checkpoint ",A)') cpfile
             end if
          end if
       else
          call sys(isys)%f(ifield)%write_chk_cps(cpfile,errmsg)
          if (len_trim(errmsg) > 0) then
             write (uout,'("!! Warning !! Could not write the checkpoint ",A,": ",A)') cpfile, errmsg
          else
             write (uout,'("* Critical points written to the checkpoint: ",A/)') cpfile
          end if
       end if
    end if

    if (cancelled) then
       w%errmsg = "The " // merge("deletion","search  ",kind == cpjob_delete)
       w%errmsg = trim(w%errmsg) // " was cancelled: the critical points are unchanged"
       return
    elseif (.not.ok) then
       w%errmsg = "AUTO did not run: the options were rejected (see the output console)"
       return
    end if

    ! an export or an estimate leaves the CP list as it was
    if (.not.changes) then
       if (kind == cpjob_export) call okfile_save_dir(w%okfile)
       return
    end if
    call sysc(isys)%post_event(lastchange_cplist)

    ! a search from one point that found nothing new
    if (kind == cpjob_add .and. sys(isys)%f(ifield)%ncp == ncp0) &
       w%errmsg = "No new critical point: the search from that point found none, or one already in the list"

    ! the objects in the view
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
    use gui_main, only: ColorHighlightScene
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_calcheight, iw_table_column,&
       iw_table_headers_row, iw_highlight_selectable, iw_atom_button, iw_cell_right,&
       iw_close_button, iw_table_sort_specs, iw_helpermark, iw_tblwrite_none
    use systemmod, only: npointprop_keywords, pointprop_keywords, pointprop_keyword_desc
    use param, only: bohrtoa
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    integer, intent(inout) :: ihover(2)
    logical, intent(inout) :: ttshown

    ! the columns of the results table: the built-in ones, then one
    ! for each point property of the system (ic_NBUILTIN + ip - 1), up
    ! to the columns an ImGui table can have (IMGUI_TABLE_MAX_COLUMNS)
    integer(c_int), parameter :: ic_del = 0, ic_cp = 1, ic_x = 2, ic_y = 3, ic_z = 4, ic_wyc = 5,&
       ic_f = 6, ic_fval = 7, ic_grad = 8, ic_gradval = 9, ic_lap = 10, ic_lapval = 11,&
       ic_l1 = 12, ic_l2 = 13, ic_l3 = 14, ic_ellip = 15, ic_ends = 16, ic_path = 17,&
       ic_d1 = 18, ic_d2 = 19, ic_ang = 20, ic_NBUILTIN = 21
    integer, parameter :: maxtablecol = 64
    integer, parameter :: maxppcol = maxtablecol - ic_NBUILTIN

    ! the tables of the window that can be written as text
    integer, parameter :: itable_results = 1

    ! what is at the end of a bond path (path_end)
    integer, parameter :: endk_none = 0, endk_atom = 1, endk_cp = 2

    ! the headers of the built-in columns (the position and Wyckoff
    ! ones depend on the system), and the short definitions of the
    ! properties among them, for the list of properties
    character(len=14), parameter :: colhdr(0:ic_NBUILTIN-1) = (/character(len=14) ::&
       "(delete)","CP","x","y","z","Wyc","Field","Field (val)","|Gradient|","|Grad| (val)",&
       "Laplacian","Lap. (val)","λ1","λ2","λ3","Ellipticity","Endpoints","Path (Å)",&
       "Dist. 1 (Å)","Dist. 2 (Å)","Angle (°)"/)
    character(len=40), parameter :: coldef(ic_f:ic_NBUILTIN-1) = (/character(len=40) ::&
       "field value","valence field value","gradient norm","valence gradient norm",&
       "Laplacian","valence Laplacian","Hessian eigenvalue 1","Hessian eigenvalue 2",&
       "Hessian eigenvalue 3","bond ellipticity (λ1/λ2 - 1)","ends of the bond path",&
       "bond path length","distance to end 1 of the bond path",&
       "distance to end 2 of the bond path","angle between the ends of the bond path"/)

    integer :: n, k, ifield, ihnuc(2), i, icp, iwrite
    integer(c_int) :: itable, flags, fixed, ncol, jc
    logical :: ch, cell, ismol, havespg, rclick
    character(kind=c_char,len=:), allocatable, target :: str1
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    type(ImVec2) :: sz

    ifield = w%cp%ifield
    if (.not.field_cps_or_say(isys,ifield)) return
    call update_table_caches(w,isys)
    call pp_sync()
    ismol = sys(isys)%c%ismolecule
    havespg = (sys(isys)%c%havesym > 0 .and. sys(isys)%c%spgavail)
    ihnuc = 0

    associate(f => sys(isys)%f(ifield))
      ! the list: the symmetry-unique CPs, or the cell CPs if there are
      ! more (the combo items are in the order of the tablecell values)
      cell = .false.
      if (f%ncpcel > f%ncp) then
         itable = int(w%cp%tablecell,c_int)
         call iw_combo_simple("Critical point list##cptablecombo","Symmetry-unique" // c_null_char //&
            "Cell" // c_null_char,itable,changed=ch,tooltips=&
            "The symmetry-unique critical points, one row each (their properties are those of&
            & all their copies)" // c_null_char //&
            "Every critical point in the unit cell, with the symmetry-unique one it is a copy of" //&
            c_null_char,ttshown=ttshown)
         if (ch) then
            w%cp%tablecell = int(itable)
            w%cp%sortdirty = .true.
         end if
         cell = (w%cp%tablecell == 1)
      end if
      n = merge(f%ncpcel,f%ncp,cell)

      ! the properties, which can be added to the table (they may
      ! change the number of columns)
      call draw_props_section()
      ncol = ic_NBUILTIN + int(min(sys(isys)%npropp,maxppcol),c_int)

      call draw_cp_summary(isys,ifield)
      call iw_helpermark("Right-click the header of the table, or open the list of properties&
         & above, to show more properties of the critical points or hide columns. Click a header&
         & to sort the table by that column.",sameline=.true.)

      ! the columns can be shown and hidden by right-clicking the
      ! header (with our own menu, which has tooltips), and the rows
      ! sorted by clicking it
      flags = ImGuiTableFlags_None
      flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
      flags = ior(flags,ImGuiTableFlags_RowBg)
      flags = ior(flags,ImGuiTableFlags_Borders)
      flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
      flags = ior(flags,ImGuiTableFlags_ScrollX)
      flags = ior(flags,ImGuiTableFlags_ScrollY)
      flags = ior(flags,ImGuiTableFlags_Sortable)
      ! imgui keeps the state of a column (sort, width) by position:
      ! a new table when a property other than the last one is
      ! removed or replaced (see pp_sync)
      str1 = "##tablecpresults" // string(w%cp%tablegen) // c_null_char
      ! the rest of the window, as the tables of the geometry window,
      ! leaving room for the selection buttons, a message, and Close
      call igGetContentRegionAvail(sz)
      sz%x = 0._c_float
      sz%y = sz%y - iw_calcheight(3,0,.true.)
      ! but at least a few rows. The window grows to fit its content
      ! (opening the list of properties makes it taller): the table's
      ! stretch beyond those rows is elastic, not content to fit
      sz%y = max(sz%y,iw_calcheight(7,0,.false.))
      w%heightslack = sz%y - iw_calcheight(7,0,.false.)
      if (igBeginTable(c_loc(str1),ncol,flags,sz,0._c_float)) then
         ! a disabled column is not drawn: the hidden ones, and those
         ! that cannot be shown
         do jc = 0, ncol-1
            fixed = ImGuiTableColumnFlags_WidthFixed
            if (jc == ic_del) fixed = ior(fixed,ior(ImGuiTableColumnFlags_NoHeaderLabel,&
               ImGuiTableColumnFlags_NoSort))
            if (jc == ic_ends) fixed = ior(fixed,ImGuiTableColumnFlags_NoSort)
            if (.not.column_on(jc)) fixed = ior(fixed,ImGuiTableColumnFlags_Disabled)
            call iw_table_column(col_name(jc) // "##cpcol" // string(jc),id=jc,flags=fixed)
         end do
         call iw_table_headers_row(freezetop=.true.,autofit=.true.,rclicked=rclick)
         iwrite = iw_tblwrite_none
         call header_menu(rclick,iwrite)

         ! the point properties in the table, evaluated at the CPs
         do jc = ic_NBUILTIN, ncol-1
            if (column_on(jc)) call pp_eval(jc-ic_NBUILTIN+1)
         end do

         ! the row order (by the CP column if none is sorted)
         if (iw_table_sort_specs(w%sortcid,w%sortdir,ic_cp,.true.)) w%cp%sortdirty = .true.
         ch = w%cp%sortdirty .or. .not.allocated(w%iord)
         if (.not.ch) ch = (size(w%iord) /= n)
         if (ch) call sort_rows()

         ! the rows; if the table is being written as text, they are
         ! captured as they are drawn (all of them)
         clipper = ImGuiListClipper_ImGuiListClipper()
         call ImGuiListClipper_Begin(clipper,n,-1._c_float)
         call w%table_write_begin(itable_results,iwrite,clipper)
         do while(ImGuiListClipper_Step(clipper))
            call c_f_pointer(clipper,clipper_f)
            do k = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
               call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
               call row_cp(w%iord(k),i,icp)
               call draw_cp_row(k,i,icp)
            end do
         end do
         call ImGuiListClipper_End(clipper)
         call ImGuiListClipper_destroy(clipper)
         call w%table_write_end(itable_results,"Critical points of field " // string(ifield) // " (" //&
            trim(f%name) // ") of system " // string(isys) // " (" // trim(sysc(isys)%seed%name) //&
            "), " // trim(merge("cell           ","symmetry-unique",cell)) // " list")
         call igEndTable()
      end if
    end associate

    ! selecting and deleting CPs
    call draw_edit_section(w,isys,iview,ttshown)

    ! a nucleus under the mouse: highlight the atom
    if (ihnuc(1) > 0) &
       call sysc(isys)%highlight_atoms(.true.,ihnuc(1:1),ihnuc(2),reshape(ColorHighlightScene,(/4,1/)))

  contains
    !> Whether column jc of the table can be shown: there are no
    !> Wyckoff positions for molecules or the cell list, valence
    !> quantities only for fields with a core contribution, and no
    !> columns for the point properties that are not a number (the
    !> stress tensor) or past the column limit.
    logical function available(jc)
      integer(c_int), intent(in) :: jc

      select case (jc)
      case (ic_wyc)
         available = .not.ismol .and. .not.cell
      case (ic_fval,ic_gradval,ic_lapval)
         available = sys(isys)%f(ifield)%usecore
      case (ic_NBUILTIN:)
         available = (jc - ic_NBUILTIN < maxppcol)
         if (available) available = (sys(isys)%propp(jc-ic_NBUILTIN+1)%ispecial == 0)
      case default
         available = .true.
      end select

    end function available

    !> The settings of column jc of the table (valid in this frame
    !> only: the window stack may be reallocated).
    function col(jc)
      integer(c_int), intent(in) :: jc
      type(cp_ppcol), pointer :: col

      if (jc < ic_NBUILTIN) then
         col => w%cp%bcol(jc)
      else
         col => w%cp%ppcol(jc-ic_NBUILTIN+1)
      end if

    end function col

    !> Whether column jc of the table is drawn.
    logical function column_on(jc)
      integer(c_int), intent(in) :: jc

      type(cp_ppcol), pointer :: c

      c => col(jc)
      column_on = c%show .and. available(jc)

    end function column_on

    !> The header of column jc of the table.
    function col_name(jc) result(str)
      integer(c_int), intent(in) :: jc
      character(len=:), allocatable :: str

      if (jc >= ic_NBUILTIN) then
         str = trim(sys(isys)%propp(jc-ic_NBUILTIN+1)%name)
      elseif (jc == ic_wyc .and. .not.havespg) then
         str = "Mul"
      elseif (ismol .and. jc >= ic_x .and. jc <= ic_z) then
         str = trim(colhdr(jc)) // " (Å)"
      else
         str = trim(colhdr(jc))
      end if

    end function col_name

    !> The description of column jc of the table, for the tooltips.
    function col_desc(jc) result(str)
      integer(c_int), intent(in) :: jc
      character(len=:), allocatable :: str

      select case (jc)
      case (ic_x,ic_y,ic_z)
         if (ismol) then
            str = "Cartesian coordinate of the critical point (Å)"
         else
            str = "Crystallographic (fractional) coordinate of the critical point"
         end if
      case (ic_wyc)
         if (havespg) then
            str = "Wyckoff position: the multiplicity (number of copies in the unit cell)&
               & and the Wyckoff letter"
         else
            str = "Multiplicity: the number of copies of the critical point in the unit cell"
         end if
      case (ic_f)
         str = "Value of the field at the critical point"
      case (ic_fval)
         str = "Value of the valence field (the field without its core contribution) at the&
            & critical point"
      case (ic_grad)
         str = "Norm of the gradient of the field at the critical point (zero up to the&
            & convergence threshold of the search)"
      case (ic_gradval)
         str = "Norm of the gradient of the valence field (the field without its core&
            & contribution) at the critical point"
      case (ic_lap)
         str = "Laplacian of the field at the critical point"
      case (ic_lapval)
         str = "Laplacian of the valence field (the field without its core contribution) at&
            & the critical point"
      case (ic_l1,ic_l2,ic_l3)
         str = "Hessian eigenvalue " // string(jc-ic_l1+1) // " of the field at the critical&
            & point (the eigenvalues in ascending order)"
      case (ic_ellip)
         str = "Bond ellipticity: the ratio of the two negative Hessian eigenvalues at the&
            & bond critical point (larger over smaller in absolute value) minus one"
      case (ic_ends)
         str = "The two nuclei or critical points at the ends of the bond path&
            & (bond critical points only)"
      case (ic_path)
         str = "Bond path length: the length of the gradient path that joins the two ends&
            & through the bond critical point (Å)"
      case (ic_d1,ic_d2)
         str = "Distance from the bond critical point to end " // string(jc-ic_d1+1) //&
            " of the bond path, in a straight line (Å)"
      case (ic_ang)
         str = "Bond angle: the angle between the directions from the bond critical point to&
            & the two ends of the bond path (°)"
      case (ic_NBUILTIN:)
         associate(pp => sys(isys)%propp(jc-ic_NBUILTIN+1))
           if (pp%ispecial == 0) then
              str = "Point property: " // trim(pp%expr)
              if (jc - ic_NBUILTIN >= maxppcol) str = str // ". The table has no room for&
                 & its column (there are too many point properties)"
           else
              str = "Point property: the stress tensor, which cannot be shown in the table"
           end if
         end associate
      case default
         str = ""
      end select

    end function col_desc

    !> The value of numeric column jc (a property) of the table for
    !> symmetry-unique CP i, in val. False if the cell is empty: no
    !> such property for this CP, or it could not be evaluated.
    logical function cell_value(i,jc,val)
      integer, intent(in) :: i
      integer(c_int), intent(in) :: jc
      real*8, intent(inout) :: val

      cell_value = .true.
      associate(f => sys(isys)%f(ifield))
        ! the bond path quantities: bond CPs with both ends found
        if (jc == ic_path .or. jc == ic_d1 .or. jc == ic_d2 .or. jc == ic_ang) then
           cell_value = f%isbcp(f%cp(i))
           if (cell_value) cell_value = all(f%cp(i)%ipath > 0)
           if (.not.cell_value) return
        end if
        select case (jc)
        case (ic_f)
           val = f%cp(i)%s%f
        case (ic_fval)
           val = f%cp(i)%s%fval
        case (ic_grad)
           val = f%cp(i)%s%gfmod
        case (ic_gradval)
           val = f%cp(i)%s%gfmodval
        case (ic_lap)
           val = f%cp(i)%s%del2f
        case (ic_lapval)
           val = f%cp(i)%s%del2fval
        case (ic_l1,ic_l2,ic_l3)
           val = f%cp(i)%s%hfeval(jc-ic_l1+1)
        case (ic_ellip)
           cell_value = f%isbcp(f%cp(i))
           if (cell_value) cell_value = (abs(f%cp(i)%s%hfeval(2)) > 0d0)
           if (cell_value) val = f%cp(i)%s%hfeval(1)/f%cp(i)%s%hfeval(2)-1d0
        case (ic_path)
           val = sum(f%cp(i)%brpathlen) * bohrtoa
        case (ic_d1,ic_d2)
           val = f%cp(i)%brdist(jc-ic_d1+1) * bohrtoa
        case (ic_ang)
           val = f%cp(i)%brang
        case (ic_NBUILTIN:)
           cell_value = (w%cp%ppstat(i,jc-ic_NBUILTIN+1) == 1)
           if (cell_value) val = w%cp%ppval(i,jc-ic_NBUILTIN+1)
        case default
           cell_value = .false.
        end select
      end associate

    end function cell_value

    !> Value val of numeric column jc as text, in the notation and
    !> with the decimal places chosen for it.
    function cell_str(jc,val) result(str)
      integer(c_int), intent(in) :: jc
      real*8, intent(in) :: val
      character(len=:), allocatable :: str

      type(cp_ppcol), pointer :: c

      c => col(jc)
      str = string(val,merge('e','f',c%expo),decimal=int(c%ndec))

    end function cell_str

    !> Keep the cache of the point property values and the column
    !> settings for the current list of point properties: if it
    !> changed (or the cache was reset), clear the values and rebuild
    !> the settings, keeping those of the properties still there (by
    !> name).
    subroutine pp_sync()
      integer :: ip, np, nold, ncp, k
      logical :: same
      type(cp_ppcol), allocatable :: newcol(:)

      np = sys(isys)%npropp
      nold = 0
      if (allocated(w%cp%ppfor)) then
         nold = size(w%cp%ppfor)
         same = (nold == np)
         if (same) same = all(w%cp%ppfor%name == sys(isys)%propp(1:np)%name .and.&
            w%cp%ppfor%expr == sys(isys)%propp(1:np)%expr)
         if (same) return
         ! not just appended: the table columns must start over
         same = (nold <= np)
         if (same) same = all(w%cp%ppfor%name == sys(isys)%propp(1:nold)%name .and.&
            w%cp%ppfor%expr == sys(isys)%propp(1:nold)%expr)
         if (.not.same) w%cp%tablegen = w%cp%tablegen + 1
      end if
      w%cp%ppfor = sys(isys)%propp(1:np)

      ! the column settings
      allocate(newcol(np))
      do ip = 1, np
         k = 0
         if (allocated(w%cp%ppcol)) k = findloc(w%cp%ppcol%name,sys(isys)%propp(ip)%name,1)
         if (k > 0) newcol(ip) = w%cp%ppcol(k)
         newcol(ip)%name = sys(isys)%propp(ip)%name
      end do
      call move_alloc(newcol,w%cp%ppcol)

      ! the values
      ncp = sys(isys)%f(ifield)%ncp
      if (allocated(w%cp%ppval)) deallocate(w%cp%ppval)
      if (allocated(w%cp%ppstat)) deallocate(w%cp%ppstat)
      allocate(w%cp%ppval(ncp,np),w%cp%ppstat(ncp,np))
      w%cp%ppstat = 0
      w%cp%sortdirty = .true.

    end subroutine pp_sync

    !> Evaluate point property ip at the symmetry-unique CPs, unless
    !> it already is.
    subroutine pp_eval(ip)
      use arithmetic, only: token, pretokenize
      use systemmod, only: system
      integer, intent(in) :: ip

      integer :: i
      character(len=:), allocatable :: errmsg
      type(token), allocatable :: toklist(:)
      type(system), pointer :: syl

      ! the nuclei are always in the list, so the first CP tells
      if (w%cp%ppstat(1,ip) /= 0) return
      w%cp%ppstat(:,ip) = -1
      if (sys(isys)%propp(ip)%ispecial /= 0) return

      ! the expression is parsed once, for all the CPs
      syl => sys(isys)
      errmsg = ""
      call pretokenize(trim(sys(isys)%propp(ip)%expr),toklist,errmsg,c_loc(syl))
      if (len_trim(errmsg) > 0) return
      associate(c => sys(isys)%c, f => sys(isys)%f(ifield))
        do i = 1, f%ncp
           errmsg = ""
           w%cp%ppval(i,ip) = sys(isys)%eval(sys(isys)%propp(ip)%expr,errmsg,c%x2c(f%cp(i)%x),toklist)
           w%cp%ppstat(i,ip) = merge(1,-1,len_trim(errmsg) == 0)
        end do
      end associate

    end subroutine pp_eval

    !> The collapsible list of the properties calculated at the
    !> critical points: the built-in ones and the point properties of
    !> the system, each with the checkbox that shows its column, and
    !> the form to add point properties (as POINTPROP).
    subroutine draw_props_section()
      use utils, only: iw_checkbox, iw_button, iw_inputtext, iw_arith_help_button,&
         iw_intstepper

      character(kind=c_char,len=:), allocatable, target :: strsec, strtab
      character(len=:), allocatable :: stropt, sttip
      integer(c_int) :: jc, flt
      integer :: ip, idel, nrow, k, ifmt
      logical :: ldum, ch, avail
      type(cp_ppcol), pointer :: c
      type(ImVec2) :: szt
      character(len=:), allocatable :: suffix

      ! a framed header across the window, so it stands out
      strsec = "Properties calculated at the critical points##cpprops" // c_null_char
      if (.not.igTreeNodeEx_Str(c_loc(strsec),ImGuiTreeNodeFlags_Framed)) return

      ! the table of properties: show, name, definition, notation and
      ! decimal places of the numbers in the table, delete
      nrow = sys(isys)%npropp
      do jc = ic_f, ic_ang
         if (available(jc)) nrow = nrow + 1
      end do
      flt = ImGuiTableFlags_None
      flt = ior(flt,ImGuiTableFlags_NoSavedSettings)
      flt = ior(flt,ImGuiTableFlags_RowBg)
      flt = ior(flt,ImGuiTableFlags_Borders)
      flt = ior(flt,ImGuiTableFlags_SizingFixedFit)
      flt = ior(flt,ImGuiTableFlags_ScrollY)
      strtab = "##tablecpprops" // c_null_char
      szt%x = 0._c_float
      szt%y = iw_calcheight(min(nrow,5)+1,0,.false.)
      idel = 0
      if (igBeginTable(c_loc(strtab),6_c_int,flt,szt,0._c_float)) then
         call iw_table_column("Show",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Property",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Definition",id=2_c_int,flags=ImGuiTableColumnFlags_WidthStretch)
         call iw_table_column("Format",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("Decimals",id=4_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
         call iw_table_column("(delete)",id=5_c_int,flags=ior(ImGuiTableColumnFlags_WidthFixed,&
            ImGuiTableColumnFlags_NoHeaderLabel))
         call iw_table_headers_row(freezetop=.true.)
         do jc = ic_f, ic_NBUILTIN+int(sys(isys)%npropp,c_int)-1
            avail = available(jc)
            if (jc < ic_NBUILTIN .and. .not.avail) cycle
            suffix = "cpprop" // string(jc)
            c => col(jc)
            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
            if (avail .and. igTableSetColumnIndex(0_c_int)) then
               ldum = iw_checkbox("##" // suffix // "show",c%show)
               call iw_tooltip("Show this property as a column of the table",ttshown)
            end if
            if (igTableSetColumnIndex(1_c_int)) then
               call iw_text(col_name(jc),alignframe=.true.)
               call iw_tooltip(col_desc(jc),ttshown)
            end if
            if (igTableSetColumnIndex(2_c_int)) then
               if (jc < ic_NBUILTIN) then
                  call iw_text(trim(coldef(jc)),alignframe=.true.)
               else
                  ip = jc - ic_NBUILTIN + 1
                  if (sys(isys)%propp(ip)%ispecial == 0) then
                     call iw_text(trim(sys(isys)%propp(ip)%expr),alignframe=.true.)
                  else
                     call iw_text("stress tensor (not a number: no column)",alignframe=.true.,&
                        disabled=.true.)
                  end if
               end if
               call iw_tooltip(col_desc(jc),ttshown)
            end if
            ! the notation and decimals of the numbers (not for the
            ! endpoints, or the properties without a column)
            avail = avail .and. jc /= ic_ends
            if (avail .and. igTableSetColumnIndex(3_c_int)) then
               ifmt = merge(1,0,c%expo)
               call iw_combo_simple("##" // suffix // "fmt","Fixed" // c_null_char // "Exp" //&
                  c_null_char,ifmt,changed=ch)
               if (ch) c%expo = (ifmt == 1)
               call iw_tooltip("Notation of the numbers of this property in the table: fixed&
                  & point (0.0123) or exponential (1.23E-02)",ttshown)
            end if
            if (avail .and. igTableSetColumnIndex(4_c_int)) then
               ldum = iw_intstepper(suffix // "dec",c%ndec,minval=0_c_int,maxval=12_c_int,&
                  ndigit=2,tooltip="Number of decimal places of this property in the table")
            end if
            if (jc >= ic_NBUILTIN) then
               if (igTableSetColumnIndex(5_c_int)) then
                  if (iw_close_button("##" // suffix // "del")) idel = jc - ic_NBUILTIN + 1
                  call iw_tooltip("Remove this point property from the system (as if it had&
                     & not been defined with POINTPROP)",ttshown)
               end if
            end if
         end do
         call igEndTable()
      end if
      if (idel > 0) then
         call sys(isys)%delete_pointprop(idel)
         call pp_sync()
      end if

      ! add a point property: an expression, or one of the keywords
      ! of POINTPROP, for the field of the window
      stropt = "Expression" // c_null_char
      sttip = "An arithmetic expression of the fields of the system" // c_null_char
      do k = 1, npointprop_keywords
         stropt = stropt // trim(pointprop_keywords(k)) // c_null_char
         sttip = sttip // trim(pointprop_keyword_desc(k)) // c_null_char
      end do
      call iw_combo_simple("##cpppkind",stropt,w%cp%ppkind,tooltips=sttip,ttshown=ttshown)
      call iw_tooltip("The kind of point property to add: an arithmetic expression, or one of&
         & the keywords of POINTPROP (calculated from the field of this window)",ttshown)
      if (w%cp%ppkind == 0) then
         call iw_text("Name",sameline=.true.)
         ldum = iw_inputtext("##cpppname",bufsize=11,textf=w%cp%ppname,width=10,sameline=.true.)
         call iw_tooltip("Name of the new property (at most 10 characters, no spaces): the&
            & header of its column",ttshown)
         call iw_text("Expression")
         ldum = iw_inputtext("##cpppexpr",bufsize=1023,textf=w%cp%ppexpr,width=36,sameline=.true.)
         call iw_tooltip("Arithmetic expression of the fields of the system (for instance,&
            & ""$1 - $2"" or ""-lag($1)""), evaluated at the critical points",ttshown)
         call iw_arith_help_button("##cpppexprhelp",ttshown)
      end if
      if (iw_button("Add##cpppadd",sameline=.true.)) call add_pointprop()
      call iw_tooltip("Add the property to the system (as POINTPROP) and show it as a column of&
         & the table. The point properties are those of the system: CPREPORT and POINT in the&
         & console report them too",ttshown)

      call igTreePop()

    end subroutine draw_props_section

    !> Add the point property of the form (w%cp%ppkind, ppname,
    !> ppexpr) to the system, and show its column. On error, the
    !> message goes to w%errmsg.
    subroutine add_pointprop()
      use tools_io, only: equali
      character(len=:), allocatable :: name, expr, errmsg
      character*10 :: name10
      integer :: k

      if (w%cp%ppkind == 0) then
         name = trim(adjustl(w%cp%ppname))
         expr = trim(adjustl(w%cp%ppexpr))
         if (len(name) == 0 .or. index(name," ") > 0) then
            w%errmsg = "The name of the property must be a single word"
            return
         end if
         if (len(expr) == 0) then
            w%errmsg = "The expression of the property is empty"
            return
         end if
         call sys(isys)%check_expression(expr,errmsg)
         if (len_trim(errmsg) > 0) then
            w%errmsg = "Invalid expression: " // trim(errmsg)
            return
         end if
      else
         ! a keyword: for the field of the window, named after the
         ! keyword (and the field, if it is not the reference)
         k = w%cp%ppkind
         name = trim(pointprop_keywords(k))
         expr = name // "(" // string(ifield) // ")"
         if (ifield /= sys(isys)%iref) name = name // "_" // string(ifield)
      end if
      if (len(name) > len(name10)) then
         w%errmsg = "The name of the property (" // name // ") is longer than 10 characters"
         return
      end if
      name10 = name
      do k = 1, sys(isys)%npropp
         if (equali(trim(sys(isys)%propp(k)%name),trim(name10))) then
            w%errmsg = "There is already a point property named " // trim(name10)
            return
         end if
      end do

      call sys(isys)%new_pointprop_string(trim(name10) // " " // expr,errmsg)
      if (len_trim(errmsg) > 0) then
         w%errmsg = "Could not add the property: " // trim(errmsg)
         return
      end if
      w%errmsg = ""
      if (w%cp%ppkind == 0) then
         w%cp%ppname = ""
         w%cp%ppexpr = ""
      end if
      call pp_sync()

    end subroutine add_pointprop

    !> The popup that shows and hides the columns, opened by
    !> right-clicking the header (rclick). It stays open while the
    !> columns are toggled. Its last entries ask (in iwrite,
    !> iw_tblwrite_*) for the table to be written as text.
    subroutine header_menu(rclick,iwrite)
      use utils, only: iw_menuitem, iw_table_write_items
      logical, intent(in) :: rclick
      integer, intent(inout) :: iwrite

      integer(c_int) :: jc
      character(kind=c_char,len=:), allocatable, target :: strpop
      type(cp_ppcol), pointer :: c

      strpop = "##cpcolumnmenu" // c_null_char
      if (rclick) call igOpenPopup_Str(c_loc(strpop),ImGuiPopupFlags_None)
      if (igBeginPopup(c_loc(strpop),ImGuiWindowFlags_None)) then
         call igPushItemFlag(ImGuiItemFlags_SelectableDontClosePopup,.true._c_bool)
         do jc = ic_x, ncol-1
            if (.not.available(jc)) cycle
            c => col(jc)
            if (iw_menuitem(col_name(jc),selected=c%show)) c%show = .not.c%show
            call iw_tooltip(col_desc(jc),ttshown)
         end do
         call igPopItemFlag()
         call igSeparator()
         call iw_table_write_items(iwrite,ttshown)
         call igEndPopup()
      end if

    end subroutine header_menu

    !> List row kr of the table (before sorting): symmetry-unique CP
    !> i, and cell CP icp in the cell list (0 otherwise).
    subroutine row_cp(kr,i,icp)
      integer, intent(in) :: kr
      integer, intent(out) :: i, icp

      i = kr
      icp = 0
      if (cell) then
         i = sys(isys)%f(ifield)%cpcel(kr)%idx
         icp = kr
      end if

    end subroutine row_cp

    !> The label of the row's CP (symmetry-unique CP i, or cell CP
    !> icp): its name, and the cell CP number in the cell list.
    function cp_label(i,icp) result(lbl)
      integer, intent(in) :: i, icp
      character(len=:), allocatable :: lbl

      lbl = trim(sys(isys)%f(ifield)%cp(i)%name)
      if (icp > 0) lbl = lbl // " " // string(icp)

    end function cp_label

    !> The Wyckoff position of symmetry-unique CP i: multiplicity and
    !> letter (only the multiplicity without a space group).
    function wyc_str(i) result(str)
      integer, intent(in) :: i
      character(len=:), allocatable :: str

      str = string(sys(isys)%f(ifield)%cp(i)%mult) // trim(wyc_letter(i))

    end function wyc_str

    !> Calculate the order of the n rows (w%iord) by column w%sortcid,
    !> in direction w%sortdir. The rows with an empty cell in that
    !> column go last in either direction, in list order. Clears the
    !> shift-click anchor, which is a position in the table.
    subroutine sort_rows()
      use tools, only: mergesort

      integer :: kr, i, icp
      integer, allocatable :: iperm(:)
      real*8, allocatable :: key(:)
      logical, allocatable :: blank(:)
      real*8 :: x(3)

      w%cp%sortdirty = .false.
      w%lastselected = 0
      allocate(key(n),blank(n))
      iperm = (/(kr, kr = 1, n)/)
      key = 0d0
      blank = .false.
      if (w%sortcid >= ic_NBUILTIN .and. w%sortcid < ncol) call pp_eval(w%sortcid-ic_NBUILTIN+1)
      associate(f => sys(isys)%f(ifield))
        do kr = 1, n
           call row_cp(kr,i,icp)
           select case (w%sortcid)
           case (ic_x,ic_y,ic_z)
              x = row_pos(i,icp,ismol)
              key(kr) = x(w%sortcid-ic_x+1)
           case (ic_wyc)
              key(kr) = 1000d0 * f%cp(i)%mult + ichar(wyc_letter(i))
           case (ic_del,ic_cp,ic_ends)
              ! the CP column: the order of the list
              key(kr) = kr
           case default
              blank(kr) = .not.cell_value(i,w%sortcid,key(kr))
           end select
        end do
      end associate
      if (n > 1) call mergesort(key,iperm,1,n)
      if (w%sortdir == 2) iperm = iperm(n:1:-1)
      w%iord = (/pack(iperm,.not.blank(iperm)),pack((/(kr, kr = 1, n)/),blank)/)

    end subroutine sort_rows

    !> Position of the row's CP (symmetry-unique CP i, or cell CP
    !> icp): Cartesian in the input frame (Å) if cart, fractional
    !> otherwise.
    function row_pos(i,icp,cart) result(x)
      integer, intent(in) :: i, icp
      logical, intent(in) :: cart
      real*8 :: x(3)

      associate(c => sys(isys)%c, f => sys(isys)%f(ifield))
        if (icp > 0) then
           x = f%cpcel(icp)%x
        else
           x = f%cp(i)%x
        end if
        if (cart) x = (c%x2c(x) + c%molx0) * bohrtoa
      end associate

    end function row_pos

    !> The Wyckoff letter of symmetry-unique CP i (blank if not known
    !> or there is no space group).
    character*1 function wyc_letter(i)
      integer, intent(in) :: i

      wyc_letter = " "
      if (.not.havespg) return
      if (i <= sys(isys)%c%nneq) then
         wyc_letter = sys(isys)%c%at(i)%wyc
      elseif (allocated(w%cp%wyc)) then
         wyc_letter = w%cp%wyc(i-sys(isys)%c%nneq)
      end if

    end function wyc_letter

    !> Row at position k of the table: symmetry-unique CP i (icp = 0)
    !> or cell CP icp, a copy of symmetry-unique CP i. The properties
    !> are those of the symmetry-unique CP.
    subroutine draw_cp_row(k,i,icp)
      integer, intent(in) :: k, i, icp

      real*8 :: x(3), xc(3), val
      character(len=:), allocatable :: suffix, lbl
      integer :: j, iu, jcp
      integer(c_int) :: jc
      logical :: clk, selrow, drawn

      associate(c => sys(isys)%c, f => sys(isys)%f(ifield))
        suffix = "_cprow" // string(merge(icp,i,icp > 0))

        ! the row selectable, which highlights the CP in the view and,
        ! clicked, selects its symmetry-unique CP for deletion (not the
        ! nuclei), and the button that deletes it (after the selectable,
        ! which lets the items after it take the clicks)
        if (igTableSetColumnIndex(ic_del)) then
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
           ! a click toggles the row and anchors a range; a shift-click
           ! selects the rows from the anchor to this one
           if (clk .and. i > c%nneq) then
              if (igIsKeyDown(ImGuiKey_ModShift) .and. w%lastselected >= 1 .and.&
                 w%lastselected <= n) then
                 do j = min(w%lastselected,k), max(w%lastselected,k)
                    call row_cp(w%iord(j),iu,jcp)
                    if (iu > c%nneq) w%cp%sel(iu) = .true.
                 end do
              else
                 w%cp%sel(i) = .not.w%cp%sel(i)
                 w%lastselected = k
              end if
           end if
           if (i > c%nneq) then
              call igSameLine(0._c_float,0._c_float)
              if (iw_close_button("##cprowdel" // suffix)) then
                 w%cp%pending_del = (/(iu == i, iu = 1, f%ncp)/)
                 call request_job(w,iview,cpjob_delete)
              end if
              if (icp > 0) then
                 call iw_tooltip("Delete this critical point, with all its copies in the cell",ttshown)
              else
                 call iw_tooltip("Delete this critical point",ttshown)
              end if
           end if
        end if

        ! the CP
        if (igTableSetColumnIndex(ic_cp)) then
           lbl = cp_label(i,icp)
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
        ! frame for molecules; the Cartesian one in the tooltip (built
        ! only when the cell is hovered)
        x = row_pos(i,icp,ismol)
        do j = 1, 3
           if (igTableSetColumnIndex(int(ic_x+j-1,c_int))) then
              call iw_cell_right(coord_str(x(j)))
              if (.not.ismol .and. igIsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByPopup)) then
                 xc = row_pos(i,icp,.true.)
                 call iw_tooltip("Cartesian: " // coord_str(xc(1)) // " " // coord_str(xc(2)) //&
                    " " // coord_str(xc(3)) // " Å",ttshown)
              end if
           end if
        end do

        ! the Wyckoff position (multiplicity and letter, as in the
        ! geometry window; only the multiplicity without a space group),
        ! with the site symmetry in the tooltip (symmetry-unique CPs of
        ! crystals only)
        if (igTableSetColumnIndex(ic_wyc)) then
           call iw_text(wyc_str(i))
           call iw_tooltip("Site symmetry: " // trim(f%cp(i)%pg),ttshown)
        end if

        ! bond CPs: the ends of the bond path
        if (f%isbcp(f%cp(i)) .and. igTableSetColumnIndex(ic_ends)) then
           call path_end_badge(i,icp,1,suffix,drawn)
           if (drawn) call igSameLine(0._c_float,-1._c_float)
           call path_end_badge(i,icp,2,suffix,drawn)
        end if

        ! the properties
        do jc = ic_f, ncol-1
           if (jc == ic_ends) cycle
           if (.not.igTableSetColumnIndex(jc)) cycle
           if (cell_value(i,jc,val)) call iw_cell_right(cell_str(jc,val))
        end do
      end associate

    end subroutine draw_cp_row

    !> The atom or CP at end j of the bond path of the BCP in row
    !> (symmetry-unique CP i, or cell CP icp): its label lbl and kind
    !> (endk_atom: atom id of list type itype; endk_cp: symmetry-unique
    !> CP id; endk_none: no CP at the end, lbl says why or is empty).
    subroutine path_end(i,icp,j,lbl,kind,id,itype)
      integer, intent(in) :: i, icp, j
      character(len=:), allocatable, intent(out) :: lbl
      integer, intent(out) :: kind, id, itype

      integer :: iend
      integer(c_int) :: idx(4)

      kind = endk_none
      id = 0
      itype = 0
      lbl = ""
      associate(c => sys(isys)%c, f => sys(isys)%f(ifield))
        if (icp > 0) then
           ! cell CP: the cell CP at the end, and its lattice vector
           iend = f%cpcel(icp)%ipath(j)
           if (iend >= 1 .and. iend <= c%ncel) then
              ! a nucleus: the cell CP is the cell atom
              idx(1) = iend
              idx(2:4) = f%cpcel(icp)%ilvec(:,j)
              kind = endk_atom
              id = iend
              itype = atlisttype_ncel_frac
              lbl = anchor_label(isys,idx,"?",species=.true.)
              return
           elseif (iend > c%ncel .and. iend <= f%ncpcel) then
              kind = endk_cp
              id = f%cpcel(iend)%idx
              lbl = trim(f%cp(id)%name) // " " // string(iend) // lvec_str(f%cpcel(icp)%ilvec(:,j))
              return
           end if
        else
           ! symmetry-unique CP: the symmetry-unique CP at the end
           iend = f%cp(i)%ipath(j)
           if (iend >= 1 .and. iend <= c%nneq) then
              kind = endk_atom
              id = iend
              itype = atlisttype_nneq
              lbl = trim(c%at(iend)%name)
              return
           elseif (iend > c%nneq .and. iend <= f%ncp) then
              kind = endk_cp
              id = iend
              lbl = trim(f%cp(iend)%name)
              return
           end if
        end if

        ! no CP at the end: the unique list tells whether the path
        ! leaves the molecule
        if (f%cp(i)%ipath(j) == -1) then
           lbl = "(leaves the molecule)"
        elseif (f%cp(i)%ipath(j) /= 0) then
           lbl = "?"
        end if
      end associate

    end subroutine path_end

    !> Badge of the atom or CP at end j of the bond path of the BCP in
    !> row (symmetry-unique CP i, or cell CP icp). drawn: whether
    !> anything was drawn.
    subroutine path_end_badge(i,icp,j,suffix,drawn)
      integer, intent(in) :: i, icp, j
      character(len=*), intent(in) :: suffix
      logical, intent(out) :: drawn

      character(len=:), allocatable :: lbl
      integer :: kind, id, itype

      call path_end(i,icp,j,lbl,kind,id,itype)
      drawn = .true.
      if (kind == endk_atom) then
         call atom_badge(id,itype,lbl // "##cpend" // string(j) // suffix)
      elseif (kind == endk_cp) then
         call cp_badge(iview,isys,ifield,id,lbl // "##cpend" // string(j) // suffix)
      elseif (len(lbl) > 0) then
         call iw_text(lbl,disabled=.true.)
      else
         drawn = .false.
      end if

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

    !> A coordinate, 4 decimals (no "-0.0000").
    function coord_str(x) result(str)
      real*8, intent(in) :: x
      character(len=:), allocatable :: str

      str = string(merge(0d0,x,abs(x) < 5d-5),'f',decimal=4)

    end function coord_str

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
  !> (CPREPORT): JSON, a VMD script, or a structure file with the CPs
  !> as extra atoms. As in Save As, the format is detected from the
  !> extension unless chosen in the combo, which then sets the
  !> extension. CPREPORT reads the format from the extension, so a
  !> chosen format that does not match it is not written.
  subroutine draw_export_tab(w,isys,iview,ttshown)
    use systems, only: sys
    use utils, only: iw_text, iw_button, iw_tooltip, iw_inputtext, iw_checkbox,&
       iw_combo_simple, iw_calcwidth, iw_dragfloat_realc, file_name_root
    use crystalmod, only: struct_detect_write_format
    use tools_io, only: lower, fopen_write, fclose, string
    use param, only: dirsep, isformat_w_unknown
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(inout) :: ttshown

    integer :: iaux, i0, i1, isformat, lu, idetect, ichosen, ifmt
    logical :: ldum, changed, ismol, userk, usecart, usecell
    character(len=:), allocatable :: ext, errexp, line

    if (.not.field_cps_or_say(isys,w%cp%ifield)) return

    ! the format combo options, on first use
    if (.not.allocated(expformat_combostr)) then
       call build_write_format_combo()
       expformat_combostr = "Auto-detect" // c_null_char // expformat_names(iexp_json) // c_null_char //&
          expformat_names(iexp_vmd) // c_null_char // write_format_combostr(icombo_fmt1:)
    end if

    ! file name: editable field plus browse button
    call iw_text("File name",highlight=.true.)
    ldum = iw_inputtext("##cpexportfile",bufsize=1023,texta=w%okfile,width=36)
    call iw_tooltip("File the critical points are written to (with auto-detect, the extension&
       & selects the format)",ttshown)
    if (iw_button("Browse...##cpexportbrowse",sameline=.true.)) &
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savecpfile,idparent=w%id,orraise=-1)
    call iw_tooltip("Choose the file with a file browser",ttshown)
    call w%okfile_warn_overwrite()

    ! format combo; if the user sets a format, change the extension to match
    call igPushItemWidth(iw_calcwidth(45,1))
    call iw_combo_simple("Format##cpexportformat",expformat_combostr,w%cp%expformat,changed=changed)
    call igPopItemWidth()
    call iw_tooltip("Format of the file: JSON, a VMD script, or a structure file with the critical&
       & points as extra atoms. If auto-detect, the format is chosen based on the file extension",ttshown)
    if (changed .and. w%cp%expformat > 0) &
       w%okfile = file_name_root(w%okfile) // "." // expformat_ext(w%cp%expformat)

    ! the format given by the file name (0 = unknown). CPREPORT needs
    ! an extension, and does not know some the detection accepts
    i0 = index(w%okfile,dirsep,back=.true.)
    i1 = index(w%okfile(i0+1:),'.',back=.true.)
    ext = ""
    if (i1 > 0) ext = lower(trim(w%okfile(i0+i1+1:)))
    if (len(ext) == 0 .or. ext == "fhi" .or. ext == "34") then
       idetect = 0
    elseif (ext == "json") then
       idetect = iexp_json
    elseif (ext == "vmd") then
       idetect = iexp_vmd
    else
       call struct_detect_write_format(trim(w%okfile),isformat)
       idetect = 0
       if (isformat /= isformat_w_unknown) idetect = iexp_fmt1 - 1 + findloc(fmtperm,isformat,1)
    end if

    ! the format written: detected or chosen
    if (w%cp%expformat == 0) then
       ichosen = idetect
       if (ichosen > 0) then
          call iw_text("Detected:",highlight=.true.)
          call iw_text(expformat_name(ichosen),sameline=.true.)
       end if
    else
       ichosen = w%cp%expformat
    end if
    ifmt = 0
    if (ichosen >= iexp_fmt1) ifmt = fmtperm(ichosen - iexp_fmt1 + 1)

    ! the options that apply to the format: the k-point grid to
    ! espresso, CASTEP and FHIaims crystals, Cartesian coordinates to
    ! FHIaims crystals, the unit cell to 3D models of crystals
    ismol = sys(isys)%c%ismolecule
    userk = (ifmt == isformat_w_qein .or. ifmt == isformat_w_castepcell .or.&
       (ifmt == isformat_w_aimsin .and. .not.ismol))
    usecart = (ifmt == isformat_w_aimsin .and. .not.ismol)
    usecell = (ifmt == isformat_w_obj .or. ifmt == isformat_w_ply .or. ifmt == isformat_w_off)&
       .and. .not.ismol
    if (ichosen > 0 .and. ichosen /= iexp_json) then
       call iw_text("Options",highlight=.true.)
       ldum = iw_checkbox("Include the gradient paths (GRAPH)##cpexpgraph",w%cp%expgraph)
       call iw_tooltip("Also write the points of the bond paths, as extra atoms",ttshown)
    end if
    if (userk) then
       ldum = iw_dragfloat_realc("RKlength##cpexprk",x1=w%cp%exprk,speed=1._c_float,&
          min=1._c_float,max=200._c_float,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Length parameter for calculating the k-point grid",ttshown)
    end if
    if (usecart) then
       ldum = iw_checkbox("Cartesian coordinates##cpexpcartesian",w%cp%expcartesian)
       call iw_tooltip("Write Cartesian atomic coordinates (atom) instead of&
          & fractional coordinates (atom_frac)",ttshown)
       call iw_text("K-point grid written to " // trim(w%okfile) // "_control (overwritten if it exists)",&
          wrap=.true.)
    end if
    if (usecell) then
       ldum = iw_checkbox("Show unit cell##cpexpdocell",w%cp%expdocell)
       call iw_tooltip("Draw the unit cell edges in the written 3D model file",ttshown)
    end if
    if (ichosen == iexp_vmd) then
       ldum = iw_checkbox("Write a pdb file (PDB)##cpexppdb",w%cp%exppdb)
       call iw_tooltip("Write the structure and the critical points to a pdb file, with the&
          & bond critical points labeled by the atoms they join, instead of an xyz file",ttshown)
       if (w%cp%exppdb) then
          ldum = iw_dragfloat_realc("Strong bond density (a.u.)##cpexpstrong",x1=w%cp%expstrong,&
             speed=0.001_c_float,min=0._c_float,max=10._c_float,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Bond critical points with a density above this value are shown as&
             & strong bonds (STRONG)",ttshown)
       end if
       call iw_text("The script reads the structure from " // file_name_root(w%okfile) //&
          merge(".pdb",".xyz",w%cp%exppdb) // ", written next to it (overwritten if it exists)",wrap=.true.)
    end if

    ! write
    errexp = ""
    if (len_trim(w%okfile) == 0) then
       errexp = "Choose a file name"
    elseif (index(trim(w%okfile)," ") > 0) then
       errexp = "The file name cannot contain spaces"
    elseif (ichosen == 0 .and. len(ext) == 0) then
       errexp = "The file name needs an extension (the format)"
    elseif (ichosen == 0) then
       errexp = "Unknown file format: ." // ext
    elseif (idetect /= ichosen) then
       errexp = "The extension does not match the selected format"
    end if
    if (iw_button("Write##cpexportwrite",danger=.true.,disabled=(len(errexp) > 0))) then
       ! the writers stop the program if the file cannot be opened:
       ! check first, opening it as they do
       lu = fopen_write(trim(w%okfile),errstop=.false.)
       if (lu < 0) then
          w%errmsg = "Cannot write the file: " // trim(w%okfile)
       else
          call fclose(lu)
          line = trim(w%okfile)
          if (ichosen /= iexp_json .and. w%cp%expgraph) line = line // " graph"
          if (userk) line = line // " " // string(real(w%cp%exprk,8),'f',decimal=4)
          if (usecart .and. w%cp%expcartesian) line = line // " cartesian"
          if (usecell .and. w%cp%expdocell) line = line // " cell"
          if (ichosen == iexp_vmd .and. w%cp%exppdb) &
             line = line // " pdb strong " // string(real(w%cp%expstrong,8),'e',decimal=6)
          call request_job(w,iview,cpjob_export,line)
       end if
    end if
    call iw_tooltip("Write the critical points to the file (the output goes to the output console)",ttshown)
    if (len(errexp) > 0) call iw_text(errexp,danger=.true.,sameline=.true.)

  end subroutine draw_export_tab

  !> Extension of entry i of the format combo of the Export tab.
  function expformat_ext(i) result(ext)
    integer, intent(in) :: i
    character(len=:), allocatable :: ext

    if (i < iexp_fmt1) then
       ext = trim(expformat_exts(i))
    else
       ext = trim(fmtext(fmtperm(i - iexp_fmt1 + 1)))
    end if

  end function expformat_ext

  !> Name of entry i of the format combo of the Export tab.
  function expformat_name(i) result(name)
    integer, intent(in) :: i
    character(len=:), allocatable :: name

    if (i < iexp_fmt1) then
       name = trim(expformat_names(i))
    else
       name = trim(fmtnames(fmtperm(i - iexp_fmt1 + 1)))
    end if

  end function expformat_name

  !> Estimate the wall time of the search from the seeds of the form
  !> (w%cp%seedx) on field ifield of system isys: run the searches from
  !> a sample of the seeds, in parallel as AUTO does, until maxsample
  !> are done or the time budget runs out, and scale the elapsed time
  !> to all the seeds. The sample follows a golden-ratio sequence, so
  !> any part of it spreads over the whole seed list (all the seeding
  !> actions). The bookkeeping of the found CPs and the tracing of the
  !> bond paths are not included. Sets w%cp%estimate. Called with sy
  !> and the reference field set (auto_enter).
  subroutine estimate_search(w,isys,ifield)
    use systems, only: sys
    use global, only: eval_next
    use tools_io, only: string
    use utils, only: duration_string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, ifield

    integer, parameter :: maxsample = 512
    real*8, parameter :: budget = 2d0 ! seconds
    real*8, parameter :: golden = 0.6180339887498949d0

    integer :: n, ns, ndone, k, i, ier, lp
    integer*8 :: c0, c1, crate
    real*8 :: gfnormeps, x(3), t, ttot
    logical :: ok

    ! the gradient norm of a CP, read as AUTO reads it
    gfnormeps = 1d-12
    if (w%cp%use_gradeps) then
       lp = 1
       ok = eval_next(gfnormeps,w%cp%gradeps,lp)
    end if

    ! time the searches from the sample
    n = size(w%cp%seedx,2)
    ns = min(n,maxsample)
    ndone = 0
    call system_clock(c0,crate)
    !$omp parallel do private(i,x,ier,c1) reduction(+:ndone) schedule(dynamic)
    do k = 1, ns
       call system_clock(c1)
       if (real(c1-c0,8) / real(crate,8) > budget) cycle
       i = 1 + int(n * modulo(k * golden,1d0))
       i = min(i,n)
       x = w%cp%seedx(:,i) - sys(isys)%c%molx0
       call sys(isys)%f(ifield)%newton(x,gfnormeps,ier)
       ndone = ndone + 1
    end do
    !$omp end parallel do
    call system_clock(c1)
    t = real(c1-c0,8) / real(crate,8)
    ttot = t * real(n,8) / real(max(ndone,1),8)

    w%cp%estimate = "Estimated time: ~" // duration_string(ttot) // " (timed " //&
       string(ndone) // " of " // string(n) // " seeds; the bond paths are not included)"

  end subroutine estimate_search

  !> Compute the seeds of the form of window w for system isys again
  !> (AUTO, with no search) if the options, the system, or its geometry
  !> changed. The computation waits for the initialization threads of
  !> the system to finish, and for the user to release the widget being
  !> edited (a drag would compute them in every frame).
  subroutine update_seeds(w,isys)
    use autocp, only: autocritic
    use systems, only: sys, sysc, sys_ready, ok_system, are_threads_running
    use global, only: iunit
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    logical :: ok, stale
    integer :: ifield
    character(len=:), allocatable :: line
    type(auto_context) :: ctx

    ! the options (empty if the form cannot make valid ones)
    line = ""
    if (size(w%cp%seed) > 0) then
       if (len(form_error(w)) == 0) line = auto_options(w,isys,iunit)
    end if

    ! the field of the window (its pseudopotential charges select the
    ! atomic grids of MESH seeds)
    ifield = sys(isys)%iref
    if (sys(isys)%goodfield(w%cp%ifield)) ifield = w%cp%ifield

    ! compute the seeds again if something changed (absolute frame);
    ! the time estimate is for the old seeds
    stale = .not.allocated(w%cp%seedline)
    if (.not.stale) stale = (w%cp%seedsys /= isys) .or. (w%cp%seedline /= line) .or.&
       (w%cp%seedfield /= ifield) .or. (w%cp%seedtime /= sysc(isys)%timelastchange_geometry)
    if (stale .and. allocated(w%cp%estimate)) deallocate(w%cp%estimate)
    if (.not.stale .or. .not.ok_system(isys,sys_ready) .or. are_threads_running() .or.&
       igIsAnyItemActive()) return
    if (allocated(w%cp%seedx)) deallocate(w%cp%seedx)
    if (len(line) > 0) then
       call auto_enter(isys,ifield,ctx)
       call autocritic(line,ok,seeds=w%cp%seedx)
       call auto_leave(isys,ctx)
       if (.not.ok .and. allocated(w%cp%seedx)) deallocate(w%cp%seedx)
       if (allocated(w%cp%seedx)) w%cp%seedx = w%cp%seedx + spread(sys(isys)%c%molx0,2,size(w%cp%seedx,2))
    end if
    w%cp%seedline = line
    w%cp%seedsys = isys
    w%cp%seedfield = ifield
    w%cp%seedtime = sysc(isys)%timelastchange_geometry
    w%cp%seednew = .true.

  end subroutine update_seeds

  !> Show the last seeds computed for system isys in the scene of view
  !> iview (the anchor of the window, which may be an alternate view):
  !> keep the shapes in the view, and build them again only if the
  !> seeds changed or the shapes are gone.
  subroutine draw_seed_preview(w,isys,iview)
    use representations, only: rep_shape, shapekind_sphere
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview

    integer :: i, n
    logical :: found
    type(rep_shape), allocatable :: shp(:)

    if (.not.allocated(w%cp%seedx) .or. w%cp%seedsys /= isys) return
    if (.not.scene_ready(iview)) return
    n = min(size(w%cp%seedx,2),maxseedshow)
    if (n == 0) return
    call win(iview)%sc%show_transient_shapes(w%id,1,found=found)
    if (found .and. .not.w%cp%seednew) return
    allocate(shp(n))
    do i = 1, n
       shp(i) = rep_shape(kind=shapekind_sphere,x1=w%cp%seedx(:,i),rad=seed_rad,rgb=seed_rgb)
    end do
    call win(iview)%sc%show_transient_shapes(w%id,1,shp)
    w%cp%seednew = .false.

  end subroutine draw_seed_preview

  !> Show the clip region of the form (CLIP) of system isys in the
  !> scene of view iview, in the style of the region of the extract
  !> window: the wireframe of the box (a parallelepiped in crystals,
  !> from the fractional corners, as AUTO reads them) or a translucent
  !> sphere.
  subroutine draw_clip_preview(w,isys,iview)
    use systems, only: sys
    use representations, only: rep_shape, shapekind_sphere, shapekind_box
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview

    integer :: i
    real*8 :: x0(3), x1(3), v(3,3)

    if (.not.scene_ready(iview)) return
    associate(c => sys(isys)%c)
      if (w%cp%iclip == 1) then
         ! the corners in order (as auto_options writes them), fractional
         x0 = c%c2x(form_to_cart(isys,min(w%cp%clipx0,w%cp%clipx1)))
         x1 = c%c2x(form_to_cart(isys,max(w%cp%clipx0,w%cp%clipx1)))
         do i = 1, 3
            v(:,i) = c%m_x2c(:,i) * (x1(i) - x0(i))
         end do
         call win(iview)%sc%show_transient_shapes(w%id,2,(/rep_shape(kind=shapekind_box,&
            x1=c%x2c(x0)+c%molx0,v=v,rad=region_edgerad,rgb=region_rgb,alpha=region_alpha)/))
      else
         x0 = form_to_cart(isys,w%cp%clipx0)
         call win(iview)%sc%show_transient_shapes(w%id,2,(/rep_shape(kind=shapekind_sphere,&
            x1=x0+c%molx0,rad=w%cp%cliprad/bohrtoa,rgb=region_rgb,alpha=region_alpha)/))
      end if
    end associate

  end subroutine draw_clip_preview

  !> Whether view iview has a scene ready to draw in.
  function scene_ready(iview) result(ok)
    integer, intent(in) :: iview
    logical :: ok

    ok = (iview >= 1 .and. iview <= nwin)
    if (ok) ok = win(iview)%isopen .and. associated(win(iview)%sc)
    if (ok) ok = (win(iview)%sc%isinit /= 0)

  end function scene_ready

  !> The single-shot search of the Search tab: a search for a
  !> critical point from a point picked in the view (the middle of a
  !> bond, or a point), with the advanced options of the Search tab.
  !> The new critical point is added to the list.
  subroutine draw_single_shot(w,iview,ttshown)
    use utils, only: iw_text, iw_button, iw_tooltip
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview
    logical, intent(inout) :: ttshown

    call iw_text("Single-Shot",highlight=.true.,alignframe=.true.)
    call iw_tooltip("Search for a critical point from one point, with the advanced options&
       & below. The new critical point is added to the list",ttshown)
    if (iw_button("Pick##cpaddpick",sameline=.true.,&
       disabled=(len(form_error(w)) > 0 .or. w%cp%picking > 0))) &
       call start_pick(w,iview,cppick_add)
    call iw_tooltip("Search for a critical point starting at a point picked in the view: the middle&
       & of a bond if one is clicked, otherwise the clicked point (on the plane through the center of&
       & the scene)",ttshown)

  end subroutine draw_single_shot

  !> Draw the field combo of the critical points window, above its
  !> tabs: the field whose critical points are searched, shown, and
  !> exported. Says so if the field is not available in system isys.
  subroutine draw_field_combo(w,isys,ttshown)
    use systems, only: sys
    use utils, only: iw_text, iw_tooltip, iw_field_combo, iw_calcwidth
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    logical, intent(inout) :: ttshown

    integer(c_int) :: ifield

    call iw_text("Field",highlight=.true.,alignframe=.true.)
    call igSameLine(0._c_float,-1._c_float)
    ifield = w%cp%ifield
    if (iw_field_combo("##cpfieldcombo",isys,ifield,width=iw_calcwidth(30,1),&
       nonestr="<field not available>")) w%cp%ifield = ifield
    call iw_tooltip("Field whose critical points are searched, shown, and exported",ttshown)
    if (.not.sys(isys)%goodfield(w%cp%ifield)) &
       call iw_text("The selected field is not available in this system",danger=.true.,wrap=.true.)

  end subroutine draw_field_combo

  !> Whether field ifield of system isys has critical points other
  !> than the nuclei; if not, say so (the Results and Export tabs).
  !> Silent if the field is not available (the field combo says so).
  function field_cps_or_say(isys,ifield) result(ok)
    use systems, only: sys
    use representations, only: field_has_cps
    use utils, only: iw_text
    integer, intent(in) :: isys, ifield
    logical :: ok

    ok = sys(isys)%goodfield(ifield)
    if (.not.ok) return
    ok = field_has_cps(isys,ifield)
    if (.not.ok) &
       call iw_text("This field has no critical points other than the nuclei (search for them in&
          & the Search tab, or load a checkpoint that has them)",disabled=.true.,wrap=.true.)

  end function field_cps_or_say

  !> The editing of the CP list, under the results table: select CPs
  !> and delete the selected ones.
  subroutine draw_edit_section(w,isys,iview,ttshown)
    use systems, only: sys
    use utils, only: iw_text, iw_button, iw_tooltip, iw_helpermark
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(inout) :: ttshown

    integer :: nsel, nnuc, iu

    ! the selection (the nuclei cannot be selected), and delete
    nnuc = sys(isys)%c%nneq
    if (iw_button("All##cpselall")) w%cp%sel(nnuc+1:) = .true.
    call iw_tooltip("Select all the critical points that are not nuclei",ttshown)
    if (iw_button("None##cpselnone",sameline=.true.)) w%cp%sel = .false.
    call iw_tooltip("Unselect all critical points",ttshown)
    if (iw_button("Toggle##cpseltoggle",sameline=.true.)) w%cp%sel(nnuc+1:) = .not.w%cp%sel(nnuc+1:)
    call iw_tooltip("Select the critical points that are not selected, and vice versa",ttshown)
    nsel = count(w%cp%sel)
    if (iw_button("Delete selected##cpdelete",sameline=.true.,disabled=(nsel == 0))) then
       w%cp%pending_del = w%cp%sel
       call request_job(w,iview,cpjob_delete)
    end if
    call iw_tooltip("Delete the selected critical points (the nuclei cannot be deleted) with all&
       & their copies in the cell, and trace the bond paths again. If only the nuclei are left,&
       & the checkpoint file of the field is removed",ttshown)
    if (iw_button("Clear##cpclearlist",danger=.true.,sameline=.true.,&
       disabled=(size(w%cp%sel) <= nnuc))) then
       ! all the critical points except the nuclei
       w%cp%pending_del = (/(iu > nnuc, iu = 1, size(w%cp%sel))/)
       call request_job(w,iview,cpjob_delete)
    end if
    call iw_tooltip("Delete all the critical points except the nuclei. This also removes&
       & the checkpoint file of the field (<field file>.chk_cps), unless writing it is&
       & disabled in the advanced options of the Search tab",ttshown)
    if (nsel > 0) call iw_text(string(nsel) // " selected",sameline=.true.)

  end subroutine draw_edit_section

  !> Start a pick in view iview of the target kind (cppick_*; iseed =
  !> the seed, for a seed position). The point to add a CP from is a
  !> bond (its midpoint) or a point; a position is an atom (its
  !> position) or a point.
  subroutine start_pick(w,iview,kind,iseed)
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview, kind
    integer, intent(in), optional :: iseed

    w%cp%picking = kind
    w%cp%pickseed = 0
    if (present(iseed)) w%cp%pickseed = iseed
    w%cp%pickview = iview
    call w%cp%pick%arm()
    if (kind == cppick_add) then
       call win(iview)%viewmode_set_forced(vm_pick_bond,"Pick a bond or a point to search for a&
          & critical point from",w%id,acceptempty=.true.)
    else
       call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick an atom or a point for the position",&
          w%id,acceptempty=.true.)
    end if

  end subroutine start_pick

  !> Forget the pending pick of window w, releasing the view.
  subroutine cancel_pick(w)
    type(window), intent(inout), target :: w

    if (w%cp%picking == 0) return
    if (w%cp%pickview >= 1 .and. w%cp%pickview <= nwin) &
       call win(w%cp%pickview)%viewmode_release_forced(w%id)
    w%cp%picking = 0

  end subroutine cancel_pick

  !> Forget a pending pick of a seed position (the seeds were deleted
  !> or reset, so its seed may be gone).
  subroutine cancel_seed_pick(w)
    type(window), intent(inout), target :: w

    if (w%cp%picking == cppick_x0 .or. w%cp%picking == cppick_x1) call cancel_pick(w)

  end subroutine cancel_seed_pick

  !> A Pick button on the line of a position widget: it starts a pick
  !> in view iview of the target kind (cppick_*; iseed = the seed). id
  !> is the ImGui id.
  subroutine pick_button(w,iview,id,kind,iseed,ttshown)
    use utils, only: iw_button, iw_tooltip
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview, kind, iseed
    character(len=*), intent(in) :: id
    logical, intent(inout) :: ttshown

    if (iw_button("Pick##" // id,sameline=.true.,disabled=(w%cp%picking > 0))) &
       call start_pick(w,iview,kind,iseed)
    call iw_tooltip("Pick the position in the view: click an atom to use its position, or empty&
       & space to use the clicked point (on the plane through the center of the scene)",ttshown)

  end subroutine pick_button

  !> Handle the pending pick commanded to view iview (system isys):
  !> the picked point (cell-frame Cartesian) goes to its target, the
  !> point to add a CP from (a picked bond gives its midpoint) or a
  !> position of the form.
  subroutine poll_pick(w,isys,iview)
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview

    integer :: istat, i
    real*8 :: xc(3)

    call view_pick_result(iview,w%id,isys,w%cp%pick,istat,xc)
    if (istat == ipick_point) then
       i = w%cp%pickseed
       select case (w%cp%picking)
       case (cppick_add)
          call request_add(w,isys,iview,xc)
       case (cppick_x0)
          if (i >= 1 .and. i <= size(w%cp%seed)) w%cp%seed(i)%x0 = form_coords(isys,xc)
       case (cppick_x1)
          if (i >= 1 .and. i <= size(w%cp%seed)) w%cp%seed(i)%x1 = form_coords(isys,xc)
       case (cppick_clipx0)
          w%cp%clipx0 = clip_coords(xc)
       case (cppick_clipx1)
          w%cp%clipx1 = clip_coords(xc)
       end select
    end if
    if (istat /= ipick_pending) w%cp%picking = 0

  contains
    !> A picked CLIP position: in crystals, wrapped to the main cell
    !> (AUTO wraps the seeds there but tests the region without
    !> periodicity); a box gets its corners sorted componentwise.
    function clip_coords(xc) result(x)
      use systems, only: sys
      real*8, intent(in) :: xc(3)
      real*8 :: x(3)

      real*8 :: xmin(3)

      x = form_coords(isys,xc)
      if (sys(isys)%c%ismolecule) return
      x = x - floor(x)
      if (w%cp%iclip == 1) then
         if (w%cp%picking == cppick_clipx0) then
            xmin = min(x,w%cp%clipx1)
            w%cp%clipx1 = max(x,w%cp%clipx1)
         else
            xmin = min(x,w%cp%clipx0)
            x = max(x,w%cp%clipx0)
            w%cp%clipx0 = xmin
            return
         end if
         x = xmin
      end if
    end function clip_coords
  end subroutine poll_pick

  !> The position xc (cell-frame Cartesian bohr) in the coordinates of
  !> the form of system isys: fractional for crystals, Cartesian Å in
  !> the frame of the input coordinates for molecules.
  function form_coords(isys,xc) result(x)
    use systems, only: sys
    integer, intent(in) :: isys
    real*8, intent(in) :: xc(3)
    real*8 :: x(3)

    associate(c => sys(isys)%c)
      if (c%ismolecule) then
         x = (xc + c%molx0) * bohrtoa
      else
         x = c%c2x(xc)
      end if
    end associate

  end function form_coords

  !> The position x in the coordinates of the form of system isys
  !> (fractional for crystals, Cartesian Å in the frame of the input
  !> coordinates for molecules) as cell-frame Cartesian bohr: the
  !> inverse of form_coords.
  function form_to_cart(isys,x) result(xc)
    use systems, only: sys
    integer, intent(in) :: isys
    real*8, intent(in) :: x(3)
    real*8 :: xc(3)

    associate(c => sys(isys)%c)
      if (c%ismolecule) then
         xc = x / bohrtoa - c%molx0
      else
         xc = c%x2c(x)
      end if
    end associate

  end function form_to_cart

  !> Request the job that adds a CP to the list by a search from xc
  !> (cell-frame Cartesian bohr), for the view iview of system isys.
  subroutine request_add(w,isys,iview,xc)
    use global, only: iunit
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    real*8, intent(in) :: xc(3)

    call request_job(w,iview,cpjob_add,auto_options(w,isys,iunit,point=form_coords(isys,xc)))

  end subroutine request_add

  !> Request a blocking job of kind kind (cpjob_*) for view iview, with
  !> the AUTO options line (search, add). The main loop runs it after
  !> the frame with its overlay; a second request in the same frame is
  !> ignored.
  subroutine request_job(w,iview,kind,line)
    type(window), intent(inout), target :: w
    integer, intent(in) :: iview, kind
    character(len=*), intent(in), optional :: line

    if (.not.w%request_block()) return
    w%cp%pending_kind = kind
    if (present(line)) w%cp%pending_line = line
    w%cp%pending_view = iview
    w%errmsg = ""

  end subroutine request_job

  !> Recompute the caches of the results table if the CP list of the
  !> field (or the field) changed: the Wyckoff letters of the
  !> symmetry-unique CPs. The selection is cleared.
  subroutine update_table_caches(w,isys)
    use systems, only: sys, sysc
    use representations, only: cp_wyckoff
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    if (w%cp%tfield == w%cp%ifield .and. w%cp%ttime == sysc(isys)%timelastchange_cplist) return
    w%cp%tfield = w%cp%ifield
    w%cp%ttime = sysc(isys)%timelastchange_cplist
    call cp_wyckoff(isys,w%cp%ifield,w%cp%wyc)
    if (allocated(w%cp%sel)) deallocate(w%cp%sel)
    allocate(w%cp%sel(sys(isys)%f(w%cp%ifield)%ncp))
    w%cp%sel = .false.
    w%cp%sortdirty = .true.
    call w%table_write_cancel()
    if (allocated(w%cp%ppfor)) deallocate(w%cp%ppfor)

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
    use representations, only: field_has_cps
    use global, only: iunit
    use gui_main, only: g
    use tools_io, only: string
    use utils, only: iw_text, iw_button, iw_tooltip, iw_checkbox, iw_helpermark
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
    logical, intent(inout) :: ttshown

    integer :: i, idel
    logical :: ldum, ok
    real(c_float) :: xcol
    character(len=:), allocatable :: xunit, errrun, str2
    character(kind=c_char,len=:), allocatable, target :: str1
    type(ImVec2) :: sz

    ! maximum height of the seed list, in rows
    integer, parameter :: maxrow_seeds = 12

    if (sys(isys)%c%ismolecule) then
       xunit = " (Å)"
    else
       xunit = " (fractional)"
    end if

    if (.not.sys(isys)%goodfield(w%cp%ifield)) return

    ! a search from one point, with the advanced options below; the new
    ! critical point is added to the list
    call draw_single_shot(w,iview,ttshown)

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
          if (draw_seed(w,isys,iview,i,xunit,xcol,ttshown)) idel = i
       end do
       w%cp%seedh = igGetCursorPosY() - g%Style%ItemSpacing%y + g%Style%WindowPadding%y
    end if
    call igEndChild()
    if (idel > 0) then
       w%cp%seed = [w%cp%seed(:idel-1), w%cp%seed(idel+1:)]
       call cancel_seed_pick(w)
    end if
    if (iw_button("Add seeds##cpseedadd")) call add_seed(w,isys,cpseed_point)
    call iw_tooltip("Add a set of seeds, the points where the CP searches are started from",ttshown)
    if (iw_button("Reset seeds##cpseedreset",sameline=.true.)) then
       call reset_seeds(w,isys)
       call cancel_seed_pick(w)
    end if
    call iw_tooltip("Go back to the seeding AUTO uses by default for this system",ttshown)
    ldum = iw_checkbox("Show seeds##cpshowseeds",w%cp%showseeds)
    call iw_tooltip("Show the starting points of the searches in the view (grey spheres), after&
       & the pruning and the clipping",ttshown)

    ! the number of seeds the search would start from
    call update_seeds(w,isys)
    if (allocated(w%cp%seedx) .and. w%cp%seedsys == isys) then
       str2 = "(" // string(size(w%cp%seedx,2)) // " seeds"
       if (w%cp%showseeds .and. size(w%cp%seedx,2) > maxseedshow) &
          str2 = str2 // ", showing " // string(maxseedshow)
       call iw_text(str2 // ")",sameline=.true.)
    end if

    ! advanced options
    call draw_advanced(w,isys,iview,xunit,ttshown)

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
    ok = (len(errrun) == 0) .and. allocated(w%cp%seedx) .and. (w%cp%seedsys == isys)
    if (ok) ok = (size(w%cp%seedx,2) > 0)
    if (iw_button("Estimate time##cpestimate",sameline=.true.,disabled=.not.ok)) &
       call request_job(w,iview,cpjob_estimate)
    call iw_tooltip("Estimate how long the search takes, by timing the searches from a sample&
       & of the seeds (at most a few seconds). The tracing of the bond paths after the search&
       & is not included",ttshown)
    if (len(errrun) > 0) then
       call iw_text(errrun,danger=.true.,sameline=.true.)
    elseif (allocated(w%cp%estimate)) then
       call iw_text(w%cp%estimate,wrap=.true.)
    end if

    ! the critical points the field has
    if (field_has_cps(isys,w%cp%ifield)) call draw_cp_summary(isys,w%cp%ifield)

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
       iw_textwidth("Points"),iw_textwidth("Position" // xunit),iw_textwidth("Mesh level"))

  end function seed_label_width

  !> Draw the widgets of seed i of window w: its kind, a button to
  !> remove it (returns true if pressed), and the parameters that kind
  !> uses, each with its label on the left and the widget at column
  !> xcol. Choosing a kind gives its parameters usable values.
  function draw_seed(w,isys,iview,i,xunit,xcol,ttshown) result(del)
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_intstepper, iw_dragfloat_real8,&
       iw_close_button
    use gui_main, only: g
    use systems, only: sys
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview, i
    character(len=*), intent(in) :: xunit
    real(c_float), intent(in) :: xcol
    logical, intent(inout) :: ttshown
    logical :: del

    character(len=:), allocatable :: suf, tt
    logical :: ldum, ch, ismol
    integer :: k, kt

    suf = "##cpseed" // string(i)
    call iw_text(string(i) // ".",alignframe=.true.)
    ! no Wigner-Seitz cell in a molecule: the list starts at the second
    ! kind (startsatone keeps typ = cpseed_*)
    ismol = sys(isys)%c%ismolecule
    k = 1
    kt = 1
    if (ismol) then
       k = index(seedkind_names,c_null_char) + 1
       kt = index(seedkind_tooltips,c_null_char) + 1
    end if
    call iw_combo_simple("##cpseedkind" // suf,seedkind_names(k:),w%cp%seed(i)%typ,sameline=.true.,&
       changed=ch,startsatone=ismol,tooltips=seedkind_tooltips(kt:),ttshown=ttshown)
    if (ch) then
       call kind_defaults(w,isys,i)
       if (w%cp%pickseed == i) call cancel_seed_pick(w)
    end if
    ! the close icon is one line of text high: centered on the combo frame
    call igSameLine(0._c_float,-1._c_float)
    call igSetCursorPosY(igGetCursorPosY() + g%Style%FramePadding%y)
    del = iw_close_button("##cpseeddel" // suf)
    call iw_tooltip("Remove this seed",ttshown)

    call igIndent(0._c_float)
    associate(s => w%cp%seed(i))
      select case(s%typ)
      case (cpseed_ws,cpseed_oh)
         tt = "Number of recursive subdivisions of the region (0 to 7; more seeds at each level)"
         call label("Subdivision level",tt)
         ldum = iw_intstepper("cpseed" // string(i) // "dep",s%depth,minval=0_c_int,maxval=int(maxdepth,c_int),tooltip=tt)
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
            ldum = iw_intstepper("cpseed" // string(i) // "nr",s%nr,minval=1_c_int,maxval=50_c_int)
         end if
      case (cpseed_sphere)
         call label("Center" // xunit)
         call x0_widget()
         call label("Radius (Å)")
         ldum = iw_dragfloat_real8(suf // "rad",x1=s%rad,speed=0.01d0,min=0.01d0,max=100d0,decimal=3)
         call label("Polar points")
         ldum = iw_intstepper("cpseed" // string(i) // "nth",s%ntheta,minval=1_c_int,maxval=100_c_int)
         call label("Azimuthal points")
         ldum = iw_intstepper("cpseed" // string(i) // "nph",s%nphi,minval=1_c_int,maxval=100_c_int)
         call label("Radial points")
         ldum = iw_intstepper("cpseed" // string(i) // "nr",s%nr,minval=1_c_int,maxval=50_c_int)
      case (cpseed_pair,cpseed_triplet)
         tt = "Only atoms closer than this distance are combined"
         call label("Maximum distance (Å)",tt)
         ldum = iw_dragfloat_real8(suf // "dist",x1=s%dist,speed=0.01d0,min=0d0,max=100d0,decimal=3)
         call iw_tooltip(tt,ttshown)
         if (s%typ == cpseed_pair) then
            tt = "Number of seeds on the segment between the two atoms"
            call label("Points per pair",tt)
            ldum = iw_intstepper("cpseed" // string(i) // "npts",s%npts,minval=1_c_int,maxval=50_c_int,tooltip=tt)
         end if
      case (cpseed_line)
         call label("Start" // xunit)
         call x0_widget()
         call label("End" // xunit)
         ldum = iw_dragfloat_real8(suf // "x1",x3=s%x1,speed=0.001d0,decimal=4)
         call pick_button(w,iview,"cppickx1" // string(i),cppick_x1,i,ttshown)
         call label("Points")
         ldum = iw_intstepper("cpseed" // string(i) // "npts",s%npts,minval=2_c_int,maxval=10000_c_int)
      case (cpseed_point)
         call label("Position" // xunit)
         call x0_widget()
      case (cpseed_mesh)
         call label("Mesh level")
         call iw_combo_simple("##cpseedmlevel" // suf,meshlevel_names,s%mlevel,startsatone=.true.,&
            tooltips=meshlevel_tooltips,ttshown=ttshown)
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
    !> The first position of the seed, with its Pick button.
    subroutine x0_widget()
      ldum = iw_dragfloat_real8(suf // "x0",x3=w%cp%seed(i)%x0,speed=0.001d0,decimal=4)
      call pick_button(w,iview,"cppickx0" // string(i),cppick_x0,i,ttshown)
    end subroutine x0_widget
  end function draw_seed

  !> The advanced options of AUTO, in a collapsed tree node; each one
  !> is used only if its checkbox is set (otherwise, AUTO's default).
  !> The labels are in one column (after the checkbox, if any), and the
  !> widgets start at another.
  subroutine draw_advanced(w,isys,iview,xunit,ttshown)
    use gui_main, only: g
    use utils, only: iw_tooltip, iw_checkbox, iw_inputtext, iw_dragfloat_real8,&
       iw_combo_simple, iw_text, iw_textwidth
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys, iview
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
          "Sphere" // c_null_char,w%cp%iclip,changed=ch,tooltips=&
          "Use all the seeds" // c_null_char //&
          "Only use the seeds inside a box, given by two opposite corners" // c_null_char //&
          "Only use the seeds inside a sphere, given by its center and radius" // c_null_char,&
          ttshown=ttshown)
       if (ch) then
          call clip_defaults(w,isys)
          if (w%cp%picking == cppick_clipx0 .or. w%cp%picking == cppick_clipx1) call cancel_pick(w)
       end if
       if (w%cp%iclip == 1) then
          call label("Box corner 1" // xunit)
          ldum = iw_dragfloat_real8("##cpclipx0",x3=w%cp%clipx0,speed=0.001d0,decimal=4)
          call pick_button(w,iview,"cppickclipx0",cppick_clipx0,0,ttshown)
          call label("Box corner 2" // xunit)
          ldum = iw_dragfloat_real8("##cpclipx1",x3=w%cp%clipx1,speed=0.001d0,decimal=4)
          call pick_button(w,iview,"cppickclipx1",cppick_clipx1,0,ttshown)
       elseif (w%cp%iclip == 2) then
          call label("Center" // xunit)
          ldum = iw_dragfloat_real8("##cpclipx0s",x3=w%cp%clipx0,speed=0.001d0,decimal=4)
          call pick_button(w,iview,"cppickclipx0s",cppick_clipx0,0,ttshown)
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
          line = line // " clip cube" // xstr("",min(w%cp%clipx0,w%cp%clipx1)) //&
             xstr("",max(w%cp%clipx0,w%cp%clipx1))
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
           case (cpseed_mesh)
              line = line // " " // trim(mesh_level_kw(s%mlevel))
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

  !> Draw the summary of the critical points of field ifield of system
  !> isys: the number of each type in the cell, and the Morse sum
  !> (crystals) or the Poincare-Hopf sum (molecules).
  subroutine draw_cp_summary(isys,ifield)
    use systems, only: sys
    use utils, only: iw_text
    use tools_io, only: string
    integer, intent(in) :: isys, ifield

    integer :: i, it, nt(0:3)

    nt = 0
    associate(f => sys(isys)%f(ifield))
      do i = 1, f%ncp
         it = f%cp(i)%typind
         if (it >= 0 .and. it <= 3) nt(it) = nt(it) + f%cp(i)%mult
      end do
    end associate
    call iw_text("Critical points",highlight=.true.)
    call iw_text("(n|b|r|c): " // string(nt(0)) // " | " // string(nt(1)) // " | " //&
       string(nt(2)) // " | " // string(nt(3)),sameline=.true.)
    if (sys(isys)%c%ismolecule) then
       call iw_text("Poincare-Hopf sum:",highlight=.true.)
    else
       call iw_text("Morse sum:",highlight=.true.)
    end if
    call iw_text(string(nt(0)-nt(1)+nt(2)-nt(3)),sameline=.true.)

  end subroutine draw_cp_summary

end submodule cp
