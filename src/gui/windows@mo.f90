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

! Routines for the molecular orbitals window.
submodule (windows) mo
  use interfaces_cimgui
  implicit none

  real(c_float), parameter :: mo_cached_rgb(3) = (/0.45_c_float,0.75_c_float,0.45_c_float/) ! cached-grid marker
  real(c_float), parameter :: mo_mark_frac = 0.35_c_float ! cached-grid mark radius, as a fraction of the row height
  real(c_float), parameter :: mo_sep_frac = 0.25_c_float ! HOMO/LUMO separator height, as a fraction of the row height
  ! The display settings, kept across a close and reopen. They are the
  ! user's choices, not properties of the system, and a window slot is
  ! recycled unpredictably, so holding them here rather than on the
  ! window is what makes them survive. Seeded from the defaults the first
  ! time any window asks
  logical :: mo_kept_valid = .false.
  integer(c_int) :: mo_kept_ieneunit = 0
  integer(c_int) :: mo_kept_ilevel = 0
  real*8 :: mo_kept_isoval = 0d0
  real(c_float) :: mo_kept_rgb(3,2) = 0._c_float
  real(c_float) :: mo_kept_alpha = 0._c_float

  character(len=2), parameter :: mo_eneunit_name(0:1) = &
     (/ character(len=2) :: "Ha", "eV" /) ! how the energy units are written
  character(len=5), parameter :: mo_spin_name(0:2) = &
     (/ character(len=5) :: "Both", "Alpha", "Beta" /) ! how a spin channel is named
  real(c_float), parameter :: mo_sep_rgb(3) = (/0.85_c_float,0.15_c_float,0.15_c_float/) ! HOMO/LUMO separator color

  ! energy-level diagram
  real(c_float), parameter :: mo_diag_occ_rgb(3) = (/0.25_c_float,0.55_c_float,0.90_c_float/) ! occupied level color
  real(c_float), parameter :: mo_diag_vir_rgb(3) = (/0.55_c_float,0.55_c_float,0.60_c_float/) ! virtual level color
  real(c_float), parameter :: mo_diag_barthick = 2.5_c_float ! level line weight (pixels)
  real(c_float), parameter :: mo_diag_selthick = 4.5_c_float ! selected level line weight (pixels)
  real(c_float), parameter :: mo_diag_bandfrac = 0.72_c_float ! fraction of a channel band the levels span
  real(c_float), parameter :: mo_diag_seppx = 4._c_float ! levels closer than this (pixels) are drawn side by side
  integer, parameter :: mo_diag_nstagmax = 6 ! most levels drawn side by side in one group
  integer, parameter :: mo_diag_nshow = 4 ! levels either side of the gap in the default view
  integer, parameter :: mo_diag_maxbars = 4000 ! most levels drawn in one frame
  integer, parameter :: mo_diag_narrmax = 8 ! most columns of electron arrows across a band
  integer, parameter :: mo_diag_arrbuf = 16 ! most levels of a bunch whose arrows are held back
  real(c_float), parameter :: mo_diag_wpref = 26._c_float ! diagram width the window opens at, in characters
  real(c_float), parameter :: mo_diag_hitpx = 6._c_float ! how close (pixels) the mouse must be to a level
  real*8, parameter :: mo_diag_occtol = 1d-6 ! occupation below which a level counts as empty
  real*8, parameter :: mo_time_budget = 2d0 ! seconds of sampling the automatic quality aims for
  integer, parameter :: mo_occ_len = 5 ! width of the occupation entry, in characters
  integer, parameter :: mo_ene_len = 10 ! width of the energy entry, in characters

contains

  !> Carry the cached orbital grids of any molecular-orbitals window
  !> over a reload of field ifield of system isys that only added the
  !> virtual orbitals. The occupied orbitals come back from the file
  !> unchanged (same coefficients, same packed indices 1..nmoocc), so
  !> their grids are still valid even though the reload counts as a
  !> change of the field set. Without this the whole cache -- possibly
  !> minutes of sampling -- would be dropped by the reload.
  module subroutine mo_cache_keep_after_reload(isys,ifield)
    use systems, only: sys, sys_init, ok_system
    integer, intent(in) :: isys
    integer, intent(in) :: ifield

    integer :: i, nmoall
    type(mo_gridcache), allocatable :: gaux(:)

    if (.not.ok_system(isys,sys_init)) return
    if (.not.sys(isys)%goodfield(ifield)) return
    nmoall = sys(isys)%f(ifield)%wfn%nmoall

    do i = 1, nwin
       if (.not.win(i)%isinit .or. win(i)%type /= wintype_mo) cycle
       if (win(i)%isys /= isys) cycle
       if (.not.allocated(win(i)%mo_cache%g)) cycle
       if (win(i)%mo_cache%ifield /= ifield) cycle

       ! the grids describe this field's occupied orbitals still
       win(i)%mo_cache%fieldgen = sys(isys)%fieldgen

       ! grow the table so the orbitals that just appeared can be cached
       if (size(win(i)%mo_cache%g,1) < nmoall) then
          allocate(gaux(nmoall))
          gaux(1:size(win(i)%mo_cache%g,1)) = win(i)%mo_cache%g
          call move_alloc(gaux,win(i)%mo_cache%g)
       end if
    end do

  end subroutine mo_cache_keep_after_reload

  !> Draw the molecular orbitals window. The window lists the MOs of
  !> the reference field of the system in its anchor view; selecting
  !> one shows it as a transient +/- isosurface pair in that view
  !> (kept alive by re-arming it every frame, and reaped when the
  !> selection or the window goes away). The Create Isosurface Object
  !> button copies the transient into a permanent representation.
  module subroutine draw_mo(w)
    use gui_main, only: g, fontsize, errmsg_linger
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sysc, sys, sys_init, ok_system, reload_field_with_virtuals
    use representations, only: representation, reptype_isosurface, repflavor_isosurface,&
       iso_isoval_mo, iso_alpha_def, iso_rgb_palette,&
       iso_grid_size, iso_region_to_box, iso_level_custom, iso_nlevel, iso_defaultlevel,&
       iso_level_label
    use wfn_private, only: wfn_rhf, wfn_uhf, wfn_rohf, wfn_spin_all, wfn_spin_alpha, wfn_spin_beta
    use types, only: id_mo_id, field_evaluation_avail, fieldeval_category_mo
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple, iw_checkbox,&
       iw_dragfloat_real8, iw_coloredit, iw_close_event, iw_setpos_bottomright,&
       iw_calcheight, iw_calcwidth, duration_string
    use tools_io, only: string
    use param, only: hartoev
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, mo_ok, goodparent, syschanged, changed, found, ldum, havecost, okgrid
    integer :: isys, iview, iref, itrep, irep, digits, iselold, nq(3), ilev
    real*8 :: tsel, tdum
    character(kind=c_char,len=:), allocatable, target :: s
    character(len=:), allocatable :: label, whynot
    type(ImVec2) :: szavail
    real(c_float) :: width, wdiag, szrow, szrowmin, szrowtab, szbut, wtab, wid(6)
    type(ImVec2) :: szdummy
    real*8 :: unitfactor, alpha8, occup, ener
    type(field_evaluation_avail) :: av

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window acts on the system shown in its anchor view
    goodparent = w%anchor(iview,isys,syschanged)

    ! initialize state
    if (w%firstpass) then
       if (.not.mo_kept_valid) then
          mo_kept_ieneunit = 0
          mo_kept_ilevel = 0 ! automatic: fit the sampling into mo_time_budget
          mo_kept_isoval = iso_isoval_mo
          mo_kept_rgb(:,1) = iso_rgb_palette(:,1) ! positive lobe: blue
          mo_kept_rgb(:,2) = iso_rgb_palette(:,4) ! negative lobe: red
          mo_kept_alpha = iso_alpha_def
          mo_kept_valid = .true.
       end if
       w%mo_ieneunit = mo_kept_ieneunit
       w%mo_ilevel = mo_kept_ilevel
       w%mo_isoval = mo_kept_isoval
       w%mo_rgb = mo_kept_rgb
       w%mo_alpha = mo_kept_alpha
    end if
    if (w%firstpass .or. syschanged) then
       w%errmsg = "" ! the error, if any, was about the previous system
       w%mo_selected = 0
       w%mo_scrolled = .false.
       w%mo_scrollto = 0
       w%mo_cache = mo_cache_state() ! the grids belong to the previous system
       w%mo_diag = mo_diagram_state() ! ditto for the level list
       w%mo_cost = mo_cost_state()
    end if

    ! initialize
    changed = .false.
    itrep = 0
    doquit = .not.goodparent
    if (.not.doquit) doquit = .not.associated(win(iview)%sc)

    ! Can the reference field of this system provide molecular orbitals?
    ! Each way of failing gets its own reason: they used to be folded
    ! together and reported as a fault of the field, which is wrong for
    ! two of the three and leaves the third with no message at all
    goodsys = ok_system(isys,sys_init)
    mo_ok = goodsys
    iref = 0
    whynot = ""
    if (.not.goodsys) then
       whynot = "No system to show. Select an initialized system in the tree"
    else
       iref = sys(isys)%iref
       mo_ok = sys(isys)%goodfield(iref)
       if (.not.mo_ok) then
          whynot = "The reference field of this system is not available"
       else
          call sys(isys)%f(iref)%eval_avail(av)
          mo_ok = av%avail(fieldeval_category_mo)
          if (.not.mo_ok) then
             whynot = "The reference field cannot provide molecular orbitals. Load a &
                &molecular wavefunction (wfn, wfx, fchk, molden) and make it the reference"
          else
             mo_ok = associated(win(iview)%sc)
             if (.not.mo_ok) whynot = "The parent view has no scene to draw the orbital in"
          end if
       end if
    end if

    ! a reloaded wavefunction (virtual orbitals read, say) is a different
    ! set of orbitals: centre the table on the new boundary again
    if (mo_ok) then
       if (w%mo_fieldgen /= sys(isys)%fieldgen) then
          w%mo_scrolled = .false.
          w%mo_fieldgen = sys(isys)%fieldgen
       end if
    end if

    ! header
    if (goodsys) then
       ! system name
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! reference field
       call iw_text("Reference field",highlight=.true.)
       if (sys(isys)%goodfield(iref)) then
          call iw_text("(" // string(iref) // ") " // trim(sys(isys)%f(iref)%name),sameline=.true.)
       end if
    end if
    if (.not.mo_ok) call iw_text(whynot,danger=.true.,wrap=.true.)

    ! Messages, dropped once they have been up long enough to read: they
    ! report an action that is over, not a state, so leaving one on the
    ! window makes it describe something that is no longer true
    if (len_trim(w%errmsg) > 0) then
       if (glfwGetTime() - w%timelast_errmsg < errmsg_linger) then
          call iw_text(w%errmsg,danger=w%mo_msgbad,wrap=.true.)
       else
          w%errmsg = ""
       end if
    end if

    if (mo_ok) then
       associate (wfn => sys(isys)%f(iref)%wfn)
         ! wavefunction description
         if (wfn%wfntyp == wfn_rhf) then
            s = "restricted"
         elseif (wfn%wfntyp == wfn_uhf) then
            s = "unrestricted"
         elseif (wfn%wfntyp == wfn_rohf) then
            s = "restricted-open"
         else
            s = "fractional occupations"
         end if
         call iw_text("Wavefunction",highlight=.true.)
         call iw_text(s // ", " // string(wfn%nmoocc) // " occupied MOs",sameline=.true.)
         if (wfn%hasvirtual) &
            call iw_text(", " // string(wfn%nmoall-wfn%nmoocc) // " virtual",sameline_nospace=.true.)
       end associate

       ! virtual orbitals: if the file has them but they were not read,
       ! offer to re-read it. The reload replaces the field, so it can
       ! not run inside an associate block that aliases it
       if (sys(isys)%f(iref)%has_unread_virtuals()) then
          if (iw_button("Reload with Virtual Orbitals",danger=.true.)) then
             call reload_field_with_virtuals(isys,iref,w%errmsg)
             if (len_trim(w%errmsg) > 0) then
                call mo_message(w,w%errmsg,.true.)
             else
                ! the occupied orbitals are unchanged: keep their grids
                call mo_cache_keep_after_reload(isys,iref)
                call mo_message(w,"Re-read the wavefunction with its virtual orbitals",.false.)
             end if
          end if
          call iw_tooltip("Re-read the wavefunction file for the reference field,&
             & this time including the virtual (unoccupied) orbitals",ttshown)
       elseif (.not.sys(isys)%f(iref)%wfn%hasvirtual) then
          call iw_text("No virtuals available")
       end if

       associate (wfn => sys(isys)%f(iref)%wfn)
         ! drop a stale selection
         if (w%mo_selected > wfn%nmoall) w%mo_selected = 0

         ! orbital table header and energy units
         call iw_text("Orbitals",highlight=.true.,alignframe=.true.)
         if (wfn%hasene) then
            call iw_combo_simple("##moeneunit","Hartree" // c_null_char // "eV" // c_null_char,&
               w%mo_ieneunit,sameline=.true.)
            call iw_tooltip("Units for the orbital energies",ttshown)
            ldum = iw_checkbox("Diagram##moshowdiag",w%mo_showdiag,sameline=.true.)
            call iw_tooltip("Show the energy-level diagram beside the orbital table",ttshown)
         end if

         if (w%mo_ieneunit == 0) then
            unitfactor = 1d0
            digits = 4
         else
            unitfactor = hartoev
            digits = 3
         end if

         ! What is selected, in words. The tables highlight the row and
         ! the diagram thickens the level, but both can be scrolled out of
         ! sight, and for an unrestricted wavefunction this is the only
         ! place that says which spin channel the orbital belongs to
         call iw_text("Selected",highlight=.true.)
         if (w%mo_selected > 0) then
            call wfn%get_mo_info(w%mo_selected,label,occup=occup,ener=ener)
            s = wfn%get_mo_name(w%mo_selected)
            if (len_trim(label) > 0) s = s // " (" // trim(label) // ")"
            s = s // ", occ. " // string(occup,'f',decimal=2)
            if (wfn%hasene) s = s // ", " // string(ener*unitfactor,'f',decimal=digits) //&
               " " // trim(mo_eneunit_name(w%mo_ieneunit))
            call iw_text(s,sameline=.true.)
         else
            call iw_text("none",disabled=.true.,sameline=.true.)
         end if

         ! The energy-level diagram, then the orbital table: one table
         ! per spin channel when the wavefunction has two of them
         ! (unrestricted), a single combined table otherwise. Every panel
         ! is a group with exactly one heading line, so they line up, and
         ! they are all given the same height.
         ! The panels take the height the window has left, with the old
         ! ten rows as the floor: hard-coding the height meant a taller
         ! window bought nothing but empty space, which is worst exactly
         ! when there are enough orbitals to want a taller window. The
         ! reserve is counted, not measured from the laid-out content --
         ! measuring it would make it depend on the thing it is sizing.
         ! Below the panels: the Display heading, its four framed rows,
         ! the Create button, and the Close row
         call igGetContentRegionAvail(szavail)
         szrowmin = iw_calcheight(10,0,.false.)
         szrow = max(szavail%y - iw_calcheight(6,1,.true.) - g%Style%ItemSpacing%y,szrowmin)
         w%heightslack = szrow - szrowmin
         if (w%mo_showdiag .and. wfn%hasene) then
            ! The tables are as wide as their own columns, so the diagram
            ! takes everything they leave rather than a fixed share: a
            ! fixed share would strand the difference as dead space to
            ! the right of a single table, at any window width. Never
            ! narrower than its own toolbar, though.
            if (wfn%get_mo_nchannels() > 1) then
               call mo_table_widths(wfn,wfn_spin_alpha,wid,wtab)
               call mo_table_widths(wfn,wfn_spin_beta,wid,width)
               wtab = wtab + width + g%Style%ItemSpacing%x
            else
               call mo_table_widths(wfn,wfn_spin_all,wid,wtab)
            end if
            wdiag = max(szavail%x - wtab - g%Style%ItemSpacing%x,iw_calcwidth(11,2))

            ! the first time round, ask for a window just wide enough for
            ! the tables and a diagram worth looking at: how much that is
            ! depends on the wavefunction, which init_window cannot know
            if (w%firstpass) &
               w%needwidth = wtab + mo_diag_wpref * fontsize%x + g%Style%ItemSpacing%x +&
                  2._c_float * g%Style%WindowPadding%x
            call igBeginGroup()
            call iw_text("Diagram",highlight=.true.)
            call mo_draw_diagram(w,wfn,isys,iref,wdiag,szrow,unitfactor,digits,changed)
            call igEndGroup()
            call igSameLine(0._c_float,-1._c_float)
            szavail%x = szavail%x - wdiag - g%Style%ItemSpacing%x
         else
            w%mo_diag%ipress = 0 ! no diagram this frame, so no press can be pending
         end if
         ! The table panel carries a button row under its heading, in the
         ! same place the diagram carries Frontier/All, so the table is
         ! shortened by that row to keep the two panels level at the
         ! bottom. A split pair shares the one button -- it re-centers
         ! both channels -- so the beta panel spaces its row instead
         szbut = igGetFrameHeight() + g%Style%ItemSpacing%y
         szrowtab = max(szrow - szbut,4._c_float*fontsize%y)
         szdummy%x = 0._c_float
         szdummy%y = igGetFrameHeight()

         iselold = w%mo_selected
         if (wfn%get_mo_nchannels() > 1) then
            ! neither table is capped: each is drawn at the width of its
            ! own columns. Capping the first at half the room would clip
            ! its energy column whenever the two channels differ in width
            call igBeginGroup()
            call iw_text(trim(mo_spin_name(wfn_spin_alpha)),highlight=.true.)
            call mo_frontier_button(w,ttshown)
            call mo_draw_table(w,wfn,wfn_spin_alpha,0._c_float,szrowtab,unitfactor,digits,changed)
            call igEndGroup()
            call igSameLine(0._c_float,-1._c_float)
            call igBeginGroup()
            call iw_text(trim(mo_spin_name(wfn_spin_beta)),highlight=.true.)
            call igDummy(szdummy)
            call mo_draw_table(w,wfn,wfn_spin_beta,0._c_float,szrowtab,unitfactor,digits,changed)
            call igEndGroup()
         else
            call igBeginGroup()
            call iw_text("Orbitals",highlight=.true.)
            call mo_frontier_button(w,ttshown)
            call mo_draw_table(w,wfn,wfn_spin_all,0._c_float,szrowtab,unitfactor,digits,changed)
            call igEndGroup()
         end if

         ! the tables have had their chance to honor a scroll-to request
         ! from the diagram; it does not outlive the frame that made it
         w%mo_scrollto = 0

         ! the table moved the selection: bring it into the diagram's view
         if (w%mo_selected /= iselold .and. w%mo_showdiag .and. wfn%hasene) &
            call mo_diagram_seek(w)

         ! display options
         call iw_text("Display",highlight=.true.)

         ! Grid quality. The automatic setting comes first and is the
         ! default: on a large wavefunction a named level can mean a
         ! minutes-long freeze, since the sampling runs on the render
         ! thread. Index 0 is automatic, 1 to iso_nlevel the named levels
         call iw_text("Quality",alignframe=.true.)

         ! Each option carries what it would cost, so that an expensive
         ! grid is turned down before it is paid for rather than after:
         ! picking a level resamples immediately, and on a large
         ! wavefunction that is a minutes-long freeze. The costs appear
         ! once the box and the per-point cost are known, which is after
         ! the first orbital -- and the first one is safe by construction,
         ! since Automatic is the default.
         havecost = w%mo_cost%secs > 0d0
         if (havecost) havecost = w%mo_cost%matches(isys,iref)
         tsel = 0d0
         s = "Automatic"
         if (havecost) then
            tsel = mo_quality_cost(w,isys,0,nq)
            s = s // " (~" // duration_string(tsel) // ")"
         end if
         s = s // c_null_char
         do ilev = 1, iso_nlevel
            s = s // iso_level_label(ilev)
            if (havecost) then
               tdum = mo_quality_cost(w,isys,ilev,nq)
               s = s // " (~" // duration_string(tdum) // ")"
               if (w%mo_ilevel == ilev) tsel = tdum
            end if
            s = s // c_null_char
         end do
         call iw_combo_simple("##molevel",s,w%mo_ilevel,sameline=.true.,changed=ldum)
         changed = changed .or. ldum
         call iw_tooltip("Density of the grid used to calculate the MO isosurface. Automatic &
            &uses the default quality, coarsened if sampling one orbital would take longer &
            &than about " // string(nint(mo_time_budget)) // " seconds",ttshown)

         ! a choice well past the budget is worth saying out loud
         if (havecost .and. tsel > 2d0 * mo_time_budget) then
            call iw_text("slow",danger=.true.,sameline=.true.)
            call iw_tooltip("Sampling one orbital at this quality takes about " //&
               duration_string(tsel) // ", and the window does not respond while it runs",ttshown)
         end if

         ! isovalue
         if (iw_dragfloat_real8("Isovalue (a.u.)##moisoval",x1=w%mo_isoval,speed=0.001d0,&
            min=0.0001d0,max=10d0,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp)) &
            changed = .true.
         call iw_tooltip("Isovalue of the positive/negative molecular orbital isosurface pair",ttshown)

         ! lobe colors and opacity
         if (iw_coloredit("+ lobe##mocolorpos",rgb=w%mo_rgb(:,1))) changed = .true.
         call iw_tooltip("Color of the positive lobe of the orbital",ttshown)
         if (iw_coloredit("- lobe##mocolorneg",rgb=w%mo_rgb(:,2),sameline=.true.)) changed = .true.
         call iw_tooltip("Color of the negative lobe of the orbital",ttshown)
         alpha8 = real(w%mo_alpha,8)
         if (iw_dragfloat_real8("Opacity##moopacity",x1=alpha8,speed=0.01d0,min=0d0,max=1d0,&
            decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)) then
            w%mo_alpha = real(alpha8,c_float)
            changed = .true.
         end if
         call iw_tooltip("Opacity of the orbital isosurfaces (1 = opaque)",ttshown)

         ! the transient isosurface representation: created or re-armed
         ! every frame an orbital is selected; when the selection (or
         ! this window) goes away it stops being armed and the scene
         ! reaps it
         if (w%mo_selected > 0) then
            call win(iview)%sc%show_transient_iso(w%id,1,itrep,found)
            if (itrep > 0) then
               associate (r => win(iview)%sc%reptrans(itrep))
                 if (.not.found .or. r%iso%ifield /= iref) then
                    ! a fresh slot, or the reference field changed: bind
                    ! the field and make the +/- isosurface pair
                    call r%iso%set_field(isys,iref)
                    r%name = "Molecular Orbital"
                    call r%iso%add_iso()
                    changed = .true.
                 end if
                 if (r%iso%imosel /= id_mo_id .or. r%iso%imoidx /= w%mo_selected) then
                    r%iso%imosel = id_mo_id
                    r%iso%imoidx = w%mo_selected
                    changed = .true.
                 end if
                 if (changed) then
                    ! push the window options onto the representation
                    r%iso%slot(1)%isoval = w%mo_isoval
                    r%iso%slot(1)%rgb = w%mo_rgb(:,1)
                    r%iso%slot(1)%alpha = w%mo_alpha
                    r%iso%slot(2)%isoval = -w%mo_isoval
                    r%iso%slot(2)%rgb = w%mo_rgb(:,2)
                    r%iso%slot(2)%alpha = w%mo_alpha
                    win(iview)%sc%forcebuildlists = .true.
                 end if

                 ! The grid, tested every frame rather than only when
                 ! something else changed: a field reload or a geometry
                 ! edit moves the box and the per-point cost under us and
                 ! sets nothing, and the quality combo prices with both.
                 ! The automatic quality resolves to a custom grid, so
                 ! r%iso%ilevel cannot stand in for the window setting,
                 ! and isgenerated stays true once anything is built
                 if (w%mo_cost%ilevel_built /= w%mo_ilevel .or. .not.r%iso%isgenerated(isys) .or.&
                    .not.w%mo_cost%matches(isys,iref)) then
                    call commit_grid(win(iview)%sc%reptrans(itrep),okgrid)
                    if (okgrid) win(iview)%sc%forcebuildlists = .true.
                 end if

                 ! keep the per-orbital grid cache in step: save what the
                 ! renderer sampled, and hand back a grid it already has
                 ! for the orbital now selected
                 call w%mo_cache%sync(isys,win(iview)%sc%reptrans(itrep),&
                    win(iview)%sc%forcebuildlists)
               end associate
            end if
         end if

         ! create a permanent isosurface object from the transient one
         if (iw_button("Create Isosurface Object",disabled=(w%mo_selected == 0 .or. itrep == 0))) then
            call win(iview)%sc%add_representation(reptype_isosurface,repflavor_isosurface,id=irep)
            if (irep > 0 .and. itrep > 0) then
               call wfn%get_mo_info(w%mo_selected,label)
               win(iview)%sc%rep(irep)%iso = win(iview)%sc%reptrans(itrep)%iso
               ! named the way the orbital is written in the tables and in
               ! the field modifiers (a<n>/b<n>, or the packed index)
               s = "MO " // wfn%get_mo_name(w%mo_selected)
               if (len_trim(label) > 0) s = s // " (" // trim(label) // ")"
               win(iview)%sc%rep(irep)%name = s
               win(iview)%sc%forcebuildlists = .true.
               ! the view does not change -- the transient was already
               ! drawing this very surface -- so without a word there is
               ! nothing at all to show the press did anything
               call mo_message(w,"Created object " // s,.false.)
            else
               call mo_message(w,"Could not create the isosurface object",.true.)
            end if
         end if
         call iw_tooltip("Copy the displayed orbital into a permanent isosurface object,&
            & editable from the Objects menu of the view. Select an orbital first",ttshown,&
            whendisabled=.true.)
       end associate
    end if ! mo_ok

    ! carry the display choices to the next window that opens
    mo_kept_ieneunit = w%mo_ieneunit
    mo_kept_ilevel = w%mo_ilevel
    mo_kept_isoval = w%mo_isoval
    mo_kept_rgb = w%mo_rgb
    mo_kept_alpha = w%mo_alpha

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(5,1)

    ! final buttons: close
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window
    if (doquit) &
       call w%end()

  contains
    !> Commit the sampling grid of representation r from its staged
    !> region and level (the counterpart of the isosurface editor's
    !> Calculate grid button).
    subroutine commit_grid(r,ok)
      type(representation), intent(inout) :: r
      logical, intent(out) :: ok

      integer :: n(3)
      real*8 :: box(3,0:3), tdum
      logical :: okbox

      ! A degenerate region leaves nothing to sample. Say so: this path
      ! is retried every frame (the guard that calls us cannot clear
      ! without a box), and the only other symptom would be a quality
      ! combo that has quietly lost all of its cost estimates
      ok = .false.
      call iso_region_to_box(isys,r%iso%iregion,r%iso%rgn_x,box,okbox)
      if (.not.okbox) then
         call mo_message(w,"Could not work out the region to sample this orbital on",.true.)
         return
      end if
      w%mo_cost%box = box

      ! the cost of one sample point, on a grid of the default quality:
      ! measured whichever quality is selected, since the combo prices
      ! every one of them and not just the automatic
      n = iso_grid_size(isys,iso_defaultlevel,box=box)
      call mo_measure_cost(w,isys,iref,r,n)

      ! the grid this quality resolves to, from the routine the combo
      ! prices with, so that the advertised time and the grid actually
      ! built can never describe different things
      tdum = mo_quality_cost(w,isys,w%mo_ilevel,n)
      if (w%mo_ilevel == 0) then
         ! the automatic grid is recorded as a custom one, so that an
         ! object created from it reads back as what was drawn
         r%iso%ilevel = iso_level_custom
         r%iso%nptscustom = n
      else
         r%iso%ilevel = w%mo_ilevel
      end if
      call r%iso%apply_grid(n,r%iso%iregion,r%iso%rgn_x)
      w%mo_cost%ilevel_built = w%mo_ilevel
      ok = .true.

    end subroutine commit_grid

  end subroutine draw_mo

  !> Keep the per-orbital sampling-grid cache c, for system isys, in
  !> step with representation r: drop every grid when the field, its
  !> data, the sampling grid, or the geometry changed under them; store
  !> the samples the renderer just took; and install the cached grid of
  !> the selected orbital, so that an orbital is evaluated once and
  !> revisiting it only re-triangulates. Sets forcebuild when a cached
  !> grid was installed and the scene has to be rebuilt.
  module subroutine mo_cache_sync(c,isys,r,forcebuild)
    use systems, only: sys, sysc
    use types, only: id_mo_id
    class(mo_cache_state), intent(inout) :: c
    integer, intent(in) :: isys
    type(representation), intent(inout) :: r
    logical, intent(inout) :: forcebuild

    integer :: i, k, ilru
    logical :: ok

    ! the cache describes one field, one generation of its data, one
    ! applied grid, and one geometry
    ok = allocated(c%g)
    if (ok) ok = (c%ifield == r%iso%ifield) .and.&
       (c%fieldgen == sys(isys)%fieldgen) .and.&
       (c%timegeom == sysc(isys)%timelastchange_geometry) .and.&
       r%iso%grid_isapplied(c%n,c%iregion,c%rgn_x)
    if (.not.ok) then
       if (allocated(c%g)) deallocate(c%g)
       c%npts = 0
       c%iuse = 0
       if (sys(isys)%goodfield(r%iso%ifield) .and. all(r%iso%nptsxyz > 0)) then
          allocate(c%g(sys(isys)%f(r%iso%ifield)%wfn%nmoall))
          c%ifield = r%iso%ifield
          c%fieldgen = sys(isys)%fieldgen
          c%timegeom = sysc(isys)%timelastchange_geometry
          c%n = r%iso%nptsxyz
          c%iregion = r%iso%iregion_ap
          c%rgn_x = r%iso%rgn_x_ap
       end if
    end if
    if (.not.allocated(c%g)) return

    ! store the samples the renderer took for the orbital it built,
    ! if they are still the ones the current state describes
    k = r%iso%imoidx_built
    if (k >= 1 .and. k <= size(c%g) .and. allocated(r%iso%ff)) then
       if (.not.allocated(c%g(k)%ff) .and. r%iso%imosel_built == id_mo_id .and.&
          r%iso%ifield_built == r%iso%ifield .and. r%iso%fieldgen_built == sys(isys)%fieldgen .and.&
          r%iso%time_built >= r%iso%timelastapply_grid .and.&
          r%iso%time_built >= sysc(isys)%timelastchange_geometry .and.&
          all(shape(r%iso%ff) == c%n)) then
          c%g(k)%ff = r%iso%ff
          c%g(k)%outdomain = r%iso%outdomain
          c%npts = c%npts + size(c%g(k)%ff,kind=8)
          c%iuse = c%iuse + 1
          c%g(k)%iuse = c%iuse
       end if
    end if

    ! install the cached grid of the selected orbital, sparing the
    ! renderer the sampling pass it would otherwise run
    k = r%iso%imoidx
    if (r%iso%imosel == id_mo_id .and. c%iscached(k)) then
       if (r%iso%imoidx_built /= k) then
          call r%iso%set_samples(isys,c%g(k)%ff,c%g(k)%outdomain)
          forcebuild = .true.
       end if
       c%iuse = c%iuse + 1
       c%g(k)%iuse = c%iuse
    end if

    ! keep the cache within its budget, dropping the grids that have
    ! gone longest without use (never the one on display)
    do while (c%npts > mo_cache_maxpts)
       ilru = 0
       do i = 1, size(c%g)
          if (.not.allocated(c%g(i)%ff) .or. i == r%iso%imoidx) cycle
          if (ilru == 0) then
             ilru = i
          elseif (c%g(i)%iuse < c%g(ilru)%iuse) then
             ilru = i
          end if
       end do
       if (ilru == 0) exit
       c%npts = c%npts - size(c%g(ilru)%ff,kind=8)
       deallocate(c%g(ilru)%ff)
    end do

  end subroutine mo_cache_sync

  !> Whether the sampling grid of orbital imo is in cache c.
  pure module function mo_cache_iscached(c,imo) result(ok)
    class(mo_cache_state), intent(in) :: c
    integer, intent(in) :: imo
    logical :: ok

    ok = .false.
    if (.not.allocated(c%g)) return
    if (imo < 1 .or. imo > size(c%g)) return
    ok = allocated(c%g(imo)%ff)

  end function mo_cache_iscached


  !> Draw one molecular-orbital table for wavefunction wfn: the whole
  !> packed orbital list (ispin = 0) or a single spin channel (1 =
  !> alpha, 2 = beta, only meaningful for an unrestricted wavefunction).
  !> width caps the table width (0 = no cap; the table is never wider
  !> than its own columns) and height is its height. The rows run in orbital order, lowest first, with a
  !> separator row at the channel's HOMO/LUMO gap; the table is centered
  !> on that boundary once, tracked by w%mo_scrolled, and on the orbital
  !> requested by w%mo_scrollto (the diagram's half of the selection
  !> sync) whenever that is set. The selection is the packed index
  !> w%mo_selected, so a split pair of tables can only ever have one
  !> orbital selected between them; changed is set when it moves.
  subroutine mo_draw_table(w,wfn,ispin,width,height,unitfactor,digits,changed)
    use wfn_private, only: molwfn, wfn_uhf, wfn_rohf
    use utils, only: iw_text, iw_tooltip, iw_table_column
    use tools_io, only: string, ioj_right, lower
    type(window), intent(inout) :: w
    type(molwfn), intent(in) :: wfn
    integer, intent(in) :: ispin
    real(c_float), intent(in) :: width
    real(c_float), intent(in) :: height
    real*8, intent(in) :: unitfactor
    integer, intent(in) :: digits
    logical, intent(inout) :: changed

    logical(c_bool) :: selected
    logical :: hasspin
    integer :: i, j, k, nocc, nrow, jsep, jscroll, jto, ischan, spin, ncol
    integer :: icid, iccache, iclabel, icspin, icocc, icene
    real(c_float) :: wid(6), wid_tot
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: str1, strl
    character(len=:), allocatable :: label
    type(ImVec2) :: sz0, szero, sz1, p0, p1
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    real*8 :: occup, ener

    logical, save :: ttshown = .false. ! tooltip flag

    szero%x = 0
    szero%y = 0

    ! the number of orbitals in this channel and where its HOMO/LUMO
    ! gap falls (ispin = 0 gives the whole packed list)
    ischan = max(ispin,1)
    nrow = wfn%get_mo_nspin(ispin,nocc)

    ! the spin and energy columns appear only when the wavefunction
    ! provides them; a split table needs no spin column, since the table
    ! the row is in already says which channel it belongs to
    hasspin = (ispin == 0 .and. (wfn%wfntyp == wfn_uhf .or. wfn%wfntyp == wfn_rohf))
    icid = 0
    iccache = 1
    iclabel = 2
    ncol = 3
    icspin = -1
    if (hasspin) then
       icspin = ncol
       ncol = ncol + 1
    end if
    icocc = ncol
    ncol = ncol + 1
    icene = -1
    if (wfn%hasene) then
       icene = ncol
       ncol = ncol + 1
    end if

    ! The table runs in orbital order, lowest first. When the channel
    ! has virtuals, a separator row marks its HOMO/LUMO gap and the
    ! one-shot centering puts that boundary in the middle of the
    ! table; without virtuals there is no gap, and the HOMO (the
    ! last row) is centered instead.
    jsep = 0
    if (nrow > nocc .and. nocc > 0) then
       jsep = nocc + 1
       nrow = nrow + 1
       jscroll = jsep
    else
       jscroll = max(nocc,1)
    end if

    ! the table is as wide as its columns, capped by the room the caller
    ! has given it
    call mo_table_widths(wfn,ispin,wid,wid_tot)
    if (width > 0._c_float) wid_tot = min(wid_tot,width)

    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1 = "##tablemolecularorbitals" // string(ispin) // c_null_char
    sz0%x = wid_tot
    sz0%y = height
    if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
       ! header setup, with the widths fixed above
       call iw_table_column("Id",id=icid,flags=ImGuiTableColumnFlags_WidthFixed,width=wid(1))
       call iw_table_column("",id=iccache,flags=ImGuiTableColumnFlags_WidthFixed,width=wid(2))
       call iw_table_column("Label",id=iclabel,flags=ImGuiTableColumnFlags_WidthFixed,width=wid(3))
       if (hasspin) &
          call iw_table_column("Spin",id=icspin,flags=ImGuiTableColumnFlags_WidthFixed,width=wid(4))
       call iw_table_column("Occ.",id=icocc,flags=ImGuiTableColumnFlags_WidthFixed,width=wid(5))
       if (wfn%hasene) &
          call iw_table_column("Energy",id=icene,flags=ImGuiTableColumnFlags_WidthFixed,width=wid(6))
       call igTableSetupScrollFreeze(0, 1) ! top row always visible

       ! draw the header
       call igTableHeadersRow()

       ! an empty channel draws the header only; the one-shot centering
       ! has nothing to anchor on, so retire it
       if (nrow == 0) w%mo_scrolled(ischan) = .true.

       ! draw the rows through a clipper (the wavefunction may have
       ! thousands of MOs); while the one-shot centering is pending,
       ! force its row in so it can anchor the scroll
       ! a scroll-to request from the diagram: find the row the orbital
       ! falls on in this channel, if it is in this channel at all
       jto = 0
       if (w%mo_scrollto > 0) then
          do j = 1, wfn%get_mo_nspin(ispin)
             if (wfn%get_mo_index(ispin,j) == w%mo_scrollto) then
                jto = j
                if (jsep > 0 .and. jto >= jsep) jto = jto + 1
                exit
             end if
          end do
       end if

       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)
       if (.not.w%mo_scrolled(ischan)) &
          call ImGuiListClipper_ForceDisplayRangeByIndices(clipper,jscroll-1,jscroll)
       if (jto > 0) then
          call ImGuiListClipper_ForceDisplayRangeByIndices(clipper,jto-1,jto)
          ! an explicit request wins over the one-shot centering: both
          ! would call igSetScrollHereY on the same frame and only the
          ! row drawn last would survive
          w%mo_scrolled(ischan) = .true.
       end if
       do while (ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do j = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
             ! the HOMO/LUMO gap: a thick line across the table
             if (j == jsep) then
                call igTableNextRow(ImGuiTableRowFlags_None,mo_sep_frac*igGetFrameHeight())
                call igTableSetBgColor(ImGuiTableBgTarget_RowBg0,&
                   igGetColorU32_Vec4(ImVec4(mo_sep_rgb(1),mo_sep_rgb(2),mo_sep_rgb(3),&
                   1._c_float)),-1_c_int)
                if (.not.w%mo_scrolled(ischan)) then
                   if (igTableSetColumnIndex(icid)) then
                      call igSetScrollHereY(0.5_c_float)
                      w%mo_scrolled(ischan) = .true.
                   end if
                end if
                cycle
             end if

             ! display row to orbital: the rows below the separator are
             ! shifted by it, then the channel index becomes a packed one
             k = j
             if (jsep > 0 .and. j > jsep) k = k - 1
             i = wfn%get_mo_index(ispin,k)
             if (i == 0) cycle

             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             call wfn%get_mo_info(i,label,spin,occup,ener)

             ! id
             if (igTableSetColumnIndex(icid)) then
                ! selectable; clicking the selected row again deselects it
                call igAlignTextToFramePadding()
                strl = "##selectmo" // string(i) // c_null_char
                flags = ImGuiSelectableFlags_SpanAllColumns
                flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
                selected = (w%mo_selected == i)
                if (igSelectable_Bool(c_loc(strl),selected,flags,szero)) then
                   if (w%mo_selected == i) then
                      w%mo_selected = 0
                   else
                      w%mo_selected = i
                   end if
                   changed = .true.
                end if

                ! center the table on the HOMO/LUMO boundary once
                if (.not.w%mo_scrolled(ischan) .and. j == jscroll) then
                   call igSetScrollHereY(0.5_c_float)
                   w%mo_scrolled(ischan) = .true.
                end if

                ! bring the orbital the diagram just selected into view
                if (j == jto) then
                   call igSetScrollHereY(0.5_c_float)
                   w%mo_scrollto = 0
                end if

                ! text: the packed index in the combined table, the
                ! number within the channel in a split one (which is
                ! what the a<n>/b<n> field modifiers take)
                if (ispin == 0) then
                   call iw_text(string(i),sameline=.true.)
                else
                   call iw_text(string(k),sameline=.true.)
                   call iw_tooltip("Orbital " // string(i) // " of the wavefunction",ttshown)
                end if
             end if

             ! whether this orbital's sampling grid is in the cache:
             ! a dot as tall as it is wide, filling the height of the
             ! row. Drawn on top of the row selectable (emitted in the
             ! id column, before this one) so it survives the
             ! selection highlight; the dummy reserves the width the
             ! column is auto-sized to
             if (igTableSetColumnIndex(iccache)) then
                sz1%x = 0.75_c_float * igGetFrameHeight()
                sz1%y = sz1%x
                call igGetCursorScreenPos(p0)
                call igDummy(sz1)
                if (w%mo_cache%iscached(i)) then
                   p1%x = p0%x + 0.5_c_float * sz1%x
                   p1%y = p0%y + 0.5_c_float * sz1%y
                   call ImDrawList_AddCircleFilled(igGetWindowDrawList(),p1,&
                      mo_mark_frac * sz1%x,igGetColorU32_Vec4(ImVec4(mo_cached_rgb(1),&
                      mo_cached_rgb(2),mo_cached_rgb(3),1._c_float)),0_c_int)
                end if
             end if

             ! label (HOMO-1, LUMO, ...)
             if (igTableSetColumnIndex(iclabel)) &
                call iw_text(label)

             ! spin channel
             if (hasspin) then
                if (igTableSetColumnIndex(icspin)) &
                   call iw_text(lower(trim(mo_spin_name(spin))))
             end if

             ! occupation
             if (igTableSetColumnIndex(icocc)) &
                call iw_text(string(occup,'f',length=mo_occ_len,decimal=2,justify=ioj_right))

             ! energy
             if (wfn%hasene) then
                if (igTableSetColumnIndex(icene)) &
                   call iw_text(string(ener*unitfactor,'f',length=mo_ene_len,decimal=digits,justify=ioj_right))
             end if
          end do ! clipper range
       end do ! clipper step
       call ImGuiListClipper_End(clipper)
       call ImGuiListClipper_destroy(clipper)
       call igEndTable()
    end if ! igBeginTable

  end subroutine mo_draw_table

  !> Draw the molecular-orbital energy-level diagram of window w for
  !> wavefunction wfn (field ifield of system isys): one band of
  !> horizontal level bars per spin channel, occupied levels carrying
  !> their electron arrows, in a plot of the given width and height.
  !> unitfactor and digits are the energy unit shared with the table.
  !> Clicking a level moves the selection w%mo_selected, the same
  !> packed index the table uses, and sets changed.
  subroutine mo_draw_diagram(w,wfn,isys,ifield,width,height,unitfactor,digits,changed)
    use gui_main, only: g, fontsize
    use wfn_private, only: molwfn
    use utils, only: iw_text, iw_button, iw_tooltip
    use tools_io, only: string
    type(window), intent(inout), target :: w
    type(molwfn), intent(in) :: wfn
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    real(c_float), intent(in) :: width
    real(c_float), intent(in) :: height
    real*8, intent(in) :: unitfactor
    integer, intent(in) :: digits
    logical, intent(inout) :: changed

    integer :: ic, i, j, k, m, lo, hi, nch, nbar, ihover, ihovlast, imo, ndec
    integer :: ia, ib, iarr, ncand, nelec, narrpos, nbuf
    integer :: bufn(mo_diag_arrbuf)
    integer(c_int) :: flags, xflags, coloc, colvir, coltext, coldim, colsel, colhov
    integer(c_int) :: colgap
    logical :: inplot, haveplot, bufok, haselec
    real(c_double) :: xmin, xmax, ymin, ymax, ylo, yhi
    real(c_double), target :: extx(2), exty(2)
    real(c_float) :: pxa, pya, pxb, pyb, hpx, wpx, bandw, halfw, xbc, x0, x1, py
    real(c_float) :: hx0, hx1, hpy, hband, xbest
    real(c_float) :: slotw, gap2, clear, dhover, dd, harr, colw, xa, xc, dbest
    real(c_float) :: yarr(mo_diag_narrmax), lastpy, pybot
    real(c_float) :: bufx(mo_diag_arrbuf), bufy(mo_diag_arrbuf)
    real*8 :: epx, pad, dy
    type(ImVec2) :: szplot, pmin, pmax, pa, pb, mpos, sztxt, p0
    type(c_ptr) :: dl
    character(kind=c_char,len=:), allocatable, target :: strid, strylab, strfmt, strext
    character(kind=c_char,len=:), allocatable, target :: strtxt
    character(len=:), allocatable :: label
    real*8 :: occup, ener

    logical, save :: ttshown = .false. ! tooltip flag

    ! the level list, rebuilt only when the wavefunction or the unit change
    call mo_diagram_build(w,wfn,isys,ifield,unitfactor)
    if (w%mo_diag%n == 0) then
       call iw_text("No orbitals to show",disabled=.true.)
       return
    end if
    nch = w%mo_diag%nch

    ! the two canned views
    if (iw_button("Frontier##modiagfit1")) then
       w%mo_diag%ifit = 1
       w%mo_scrolled = .false. ! one control, both panels
    end if
    call iw_tooltip("Show the region around the HOMO/LUMO gap, in the diagram and the table",&
       ttshown)
    if (iw_button("All##modiagfit2",sameline=.true.)) w%mo_diag%ifit = 2
    call iw_tooltip("Show every orbital level (double-clicking the diagram does the same)",ttshown)

    strid = "##modiagram" // c_null_char
    strext = "##modiagextent" // c_null_char
    if (w%mo_ieneunit == 0) then
       strylab = "Energy (Ha)" // c_null_char
    else
       strylab = "Energy (eV)" // c_null_char
    end if
    ! tick labels no more precise than the range on screen deserves,
    ! taken from the range the previous frame settled on
    dy = w%mo_diag%ylim(2) - w%mo_diag%ylim(1)
    if (dy > 0d0) then
       ndec = min(max(ceiling(-log10(dy)) + 2,1),8)
    else
       ndec = max(digits-1,1)
    end if
    strfmt = "%." // string(ndec) // "f" // c_null_char

    ! the plot fills what the toolbar, and the row of band names when
    ! there is more than one channel, leave of the table's height
    hband = 0._c_float
    if (nch > 1) hband = fontsize%y + g%Style%ItemSpacing%y
    szplot%x = width
    szplot%y = max(height - igGetFrameHeight() - g%Style%ItemSpacing%y - hband,&
       4._c_float*fontsize%y)
    flags = ior(ImPlotFlags_NoTitle,ImPlotFlags_NoLegend)
    flags = ior(flags,ImPlotFlags_NoMouseText)
    ihover = 0
    haveplot = .false.
    if (ipBeginPlot(c_loc(strid),szplot,flags)) then
       ! x is a layout axis only, one band per spin channel: locking it
       ! keeps the bands the same width at any zoom and leaves the wheel
       ! and the drag acting on the energy axis alone
       xflags = ior(ImPlotAxisFlags_NoDecorations,ImPlotAxisFlags_Lock)
       xflags = ior(xflags,ImPlotAxisFlags_NoMenus)
       xflags = ior(xflags,ImPlotAxisFlags_NoHighlight)
       call ipSetupAxes(c_null_ptr,c_loc(strylab),xflags,ImPlotAxisFlags_NoInitialFit)
       call ipSetupAxisLimits(ImAxis_X1,0.5d0,real(nch,8)+0.5d0,ImPlotCond_Always)
       call ipSetupAxisFormat(ImAxis_Y1,c_loc(strfmt))
       if (w%mo_diag%ifit /= 0) then
          call mo_diagram_range(w,ylo,yhi)
          call ipSetupAxisLimits(ImAxis_Y1,ylo,yhi,ImPlotCond_Always)
          w%mo_diag%ifit = 0
       end if

       ! an invisible item spanning every level, so that ImPlot's own
       ! double-click gesture fits the whole spectrum
       extx(1) = 0.5d0
       extx(2) = real(nch,8) + 0.5d0
       exty(1) = w%mo_diag%erange(1)
       exty(2) = w%mo_diag%erange(2)
       call ipSetNextLineStyle(ImVec4(0._c_float,0._c_float,0._c_float,0._c_float),1._c_float)
       call ipPlotLine(c_loc(strext),c_loc(extx),c_loc(exty),2_c_int,ImPlotLineFlags_None,0_c_int)

       ! the plot rectangle, read off the two limit corners (the setup
       ! calls above have to come first: these lock it)
       call ipGetPlotCurrentLimits(xmin,xmax,ymin,ymax)
       call ipPlotToPixels(xmin,ymin,ImAxis_X1,ImAxis_Y1,pxa,pya)
       call ipPlotToPixels(xmax,ymax,ImAxis_X1,ImAxis_Y1,pxb,pyb)
       pmin%x = min(pxa,pxb)
       pmin%y = min(pya,pyb)
       pmax%x = max(pxa,pxb)
       pmax%y = max(pya,pyb)
       hpx = pmax%y - pmin%y
       wpx = pmax%x - pmin%x
       w%mo_diag%ylim(1) = ymin
       w%mo_diag%ylim(2) = ymax

       haveplot = (hpx > 4._c_float .and. wpx > 4._c_float .and. ymax > ymin)
       if (haveplot) then
          ! the axis is linear, so the energy-to-pixel map is affine and
          ! there is no need to go through ImPlot for every level
          epx = real(hpx,8) / (ymax - ymin)
          bandw = wpx / real(nch,c_float)
          halfw = 0.5_c_float * mo_diag_bandfrac * bandw
          coloc = igGetColorU32_Vec4(ImVec4(mo_diag_occ_rgb(1),mo_diag_occ_rgb(2),&
             mo_diag_occ_rgb(3),1._c_float))
          colvir = igGetColorU32_Vec4(ImVec4(mo_diag_vir_rgb(1),mo_diag_vir_rgb(2),&
             mo_diag_vir_rgb(3),1._c_float))
          coltext = igGetColorU32_Col(ImGuiCol_Text,1._c_float)
          coldim = igGetColorU32_Col(ImGuiCol_Text,0.6_c_float)
          colsel = igGetColorU32_Col(ImGuiCol_HeaderActive,1._c_float)
          colhov = igGetColorU32_Col(ImGuiCol_ButtonHovered,1._c_float)
          colgap = igGetColorU32_Vec4(ImVec4(mo_sep_rgb(1),mo_sep_rgb(2),mo_sep_rgb(3),&
             0.45_c_float))

          harr = 0.85_c_float * fontsize%y
          ! The mouse counts as being on the diagram only when the plot
          ! is the hovered window: a bare rectangle test would also fire
          ! through any window, popup or tooltip drawn over it, and a
          ! click there would start a multi-second grid sampling.
          ! Hovering has to survive an active item, or the press that
          ! ImPlot latches for its own pan would cancel every selection.
          call igGetMousePos(mpos)
          inplot = (mpos%x >= pmin%x .and. mpos%x <= pmax%x .and. &
                    mpos%y >= pmin%y .and. mpos%y <= pmax%y)
          if (inplot) inplot = igIsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem)
          dhover = huge(1._c_float)
          ihovlast = 0
          nbar = 0

          dl = igGetWindowDrawList()
          call ImDrawList_PushClipRect(dl,pmin,pmax,.true._c_bool)
          do ic = 1, nch
             xbc = pmin%x + (real(ic,c_float) - 0.5_c_float) * bandw
             lo = w%mo_diag%ich0(ic)
             hi = w%mo_diag%ich0(ic+1) - 1
             if (hi < lo) cycle

             ! The band is divided into columns for the electron arrows,
             ! so that levels too close to stack their arrows vertically
             ! can step sideways instead of going without: yarr remembers
             ! how far down each column the last arrow was drawn. A level
             ! only ever puts its arrow somewhere on its own bar, and the
             ! columns just say whether that place is already taken. The
             ! sentinel is finite because the build traps FP overflow.
             narrpos = max(min(int(2._c_float*halfw / (1.6_c_float*fontsize%x)),mo_diag_narrmax),1)
             colw = 2._c_float * halfw / real(narrpos,c_float)
             yarr = -1.e20_c_float
             nbuf = 0
             bufok = .true.
             lastpy = -1.e20_c_float

             ! the HOMO/LUMO gap of this channel, in the color the table
             ! uses for its own separator row
             if (w%mo_diag%homo(ic) < huge(1d0) .and. w%mo_diag%lumo(ic) < huge(1d0)) then
                py = pmax%y - real((0.5d0*(w%mo_diag%homo(ic)+w%mo_diag%lumo(ic))-ymin)*epx,c_float)
                pa%x = xbc - halfw
                pa%y = py
                pb%x = xbc + halfw
                pb%y = py
                call ImDrawList_AddLine(dl,pa,pb,colgap,1.5_c_float)
             end if

             ! Runs are formed over the whole channel, not just over the
             ! part on screen: a run counted from the first level inside
             ! the view would lose the members sitting in the sliver just
             ! beyond the edge, and its size -- and with it the stagger
             ! and the merged-group count -- would change as the plot is
             ! panned. Only the drawing is restricted to the view.
             i = lo
             do while (i <= hi)
                ! the run of levels this one cannot be told apart from at
                ! the current zoom
                j = i
                do while (j < hi)
                   if ((w%mo_diag%ene(j+1)-w%mo_diag%ene(i)) * epx >= real(mo_diag_seppx,8)) exit
                   j = j + 1
                end do
                m = j - i + 1
                if (w%mo_diag%ene(j) < ymin .or. w%mo_diag%ene(i) > ymax) then
                   i = j + 1
                   cycle
                end if

                ! how much room the group has to itself, for the arrows
                clear = huge(1._c_float)
                if (i > lo) clear = min(clear,real((w%mo_diag%ene(i)-w%mo_diag%ene(i-1))*epx,c_float))
                if (j < hi) clear = min(clear,real((w%mo_diag%ene(j+1)-w%mo_diag%ene(j))*epx,c_float))

                ! a run too long to draw side by side, or one that would
                ! run past the budget, collapses into a single bar -- but
                ! never a run of one, which has to stay selectable
                if (m > mo_diag_nstagmax .or. (nbar + m > mo_diag_maxbars .and. m > 1)) then
                   py = pmax%y - real((0.5d0*(w%mo_diag%ene(i)+w%mo_diag%ene(j))-ymin)*epx,c_float)
                   x0 = xbc - halfw
                   x1 = xbc + halfw
                   call ImDrawList_AddLine(dl,ImVec2(x0,py),ImVec2(x1,py),coldim,&
                      mo_diag_barthick+1.5_c_float)
                   ! the count, but only where the group has a line of
                   ! its own: in a dense manifold these would pile up
                   ! into an unreadable column, and the tooltip says the
                   ! same thing on demand
                   if (x1 - x0 > 3._c_float * fontsize%x .and. clear >= fontsize%y) then
                      strtxt = "x" // string(m) // c_null_char
                      pa%x = x0
                      pa%y = py - fontsize%y - 1._c_float
                      call ImDrawList_AddText_Vec2(dl,pa,coldim,c_loc(strtxt),c_null_ptr)
                   end if
                   nbar = nbar + 1
                   if (inplot .and. mpos%x >= x0-2._c_float .and. mpos%x <= x1+2._c_float) then
                      dd = abs(mpos%y - py)
                      if (dd < dhover .and. dd <= mo_diag_hitpx) then
                         dhover = dd
                         ihover = i
                         ihovlast = j
                         hx0 = x0
                         hx1 = x1
                         hpy = py
                      end if
                   end if

                   ! A merged pile stands for levels whose arrows cannot
                   ! be drawn at all. When it holds occupied levels it
                   ! joins the bunch and fails it, so that a lone arrow
                   ! never appears beside a crowded pile that has none.
                   haselec = .false.
                   do k = i, j
                      if (abs(w%mo_diag%occ(k) - 2d0) < mo_diag_occtol .or.&
                          abs(w%mo_diag%occ(k) - 1d0) < mo_diag_occtol) then
                         haselec = .true.
                         exit
                      end if
                   end do
                   if (haselec) then
                      pybot = pmax%y - real((w%mo_diag%ene(i)-ymin)*epx,c_float)
                      if (abs(pybot - lastpy) >= 1.1_c_float*harr) then
                         call mo_diag_flusharrows(dl,nbuf,bufok,bufx,bufy,bufn,ic == 1,harr,&
                            coloc,0.3_c_float*fontsize%x)
                         nbuf = 0
                         yarr = -1.e20_c_float
                      end if
                      bufok = .false.
                      lastpy = pmax%y - real((w%mo_diag%ene(j)-ymin)*epx,c_float)
                   end if
                   i = j + 1
                   cycle
                end if

                ! the group drawn side by side, filling the band
                slotw = 2._c_float * halfw / real(m,c_float)
                gap2 = 0._c_float
                if (m > 1) gap2 = min(1.5_c_float,0.12_c_float*slotw)
                do k = i, j
                   x0 = xbc - halfw + real(k-i,c_float) * slotw + gap2
                   x1 = x0 + slotw - 2._c_float * gap2
                   py = pmax%y - real((w%mo_diag%ene(k)-ymin)*epx,c_float)

                   ! the bar: occupied or empty
                   if (w%mo_diag%occ(k) > mo_diag_occtol) then
                      call ImDrawList_AddLine(dl,ImVec2(x0,py),ImVec2(x1,py),coloc,mo_diag_barthick)
                   else
                      call ImDrawList_AddLine(dl,ImVec2(x0,py),ImVec2(x1,py),colvir,&
                         max(mo_diag_barthick-1._c_float,1._c_float))
                   end if

                   ! the selected orbital -- the one the view is showing
                   ! -- highlighted as the table row is
                   if (w%mo_selected == w%mo_diag%imo(k)) then
                      call ImDrawList_AddLine(dl,ImVec2(x0,py),ImVec2(x1,py),coltext,mo_diag_selthick)
                      call mo_diag_markrect(dl,x0,x1,py,colsel,1.5_c_float)
                   end if

                   ! The electrons: two arrows on a doubly occupied
                   ! level, one on a singly occupied one (up in the
                   ! alpha band, down in the beta one), and none on a
                   ! fractional occupation, where only the number means
                   ! anything. They go in the free column nearest the
                   ! middle of the bar -- a column being free when the
                   ! last arrow drawn in it is more than an arrow height
                   ! away -- so a crowded region shows its occupations
                   ! side by side instead of dropping them.
                   nelec = 0
                   if (abs(w%mo_diag%occ(k) - 2d0) < mo_diag_occtol) then
                      nelec = 2
                   elseif (abs(w%mo_diag%occ(k) - 1d0) < mo_diag_occtol) then
                      nelec = 1
                   end if
                   if (nelec > 0) then
                      ! a bunch is a chain of levels each of which would
                      ! draw its arrows over the one before; the previous
                      ! one ends here if this level clears it
                      if (abs(py - lastpy) >= 1.1_c_float*harr) then
                         call mo_diag_flusharrows(dl,nbuf,bufok,bufx,bufy,bufn,ic == 1,harr,&
                            coloc,0.3_c_float*fontsize%x)
                         nbuf = 0
                         bufok = .true.
                         yarr = -1.e20_c_float
                      end if
                      lastpy = py

                      ! the places on this bar an arrow could go, keeping
                      ! clear of the cache dot, and the band column each
                      ! of them falls in
                      if (bufok) then
                         xc = 0.5_c_float * (x0 + x1)
                         ncand = max(int((x1 - x0) / (1.4_c_float*fontsize%x)),1)
                         iarr = 0
                         xbest = xc
                         dbest = huge(1._c_float)
                         do ia = 1, ncand
                            xa = x0 + (real(ia,c_float) - 0.5_c_float) * (x1 - x0) / real(ncand,c_float)
                            ib = min(max(int((xa - xbc + halfw) / colw) + 1,1),narrpos)
                            if (abs(py - yarr(ib)) < 1.1_c_float*harr) cycle
                            if (abs(xa - xc) >= dbest) cycle
                            dbest = abs(xa - xc)
                            xbest = xa
                            iarr = ib
                         end do
                         if (iarr > 0 .and. nbuf < mo_diag_arrbuf) then
                            yarr(iarr) = py
                            nbuf = nbuf + 1
                            bufx(nbuf) = xbest
                            bufy(nbuf) = py
                            bufn(nbuf) = nelec
                         else
                            ! nowhere left to put this one: the whole
                            ! bunch goes without
                            bufok = .false.
                         end if
                      end if
                   end if

                   nbar = nbar + 1
                   if (inplot .and. mpos%x >= x0-2._c_float .and. mpos%x <= x1+2._c_float) then
                      dd = abs(mpos%y - py)
                      if (dd < dhover .and. dd <= mo_diag_hitpx) then
                         dhover = dd
                         ihover = k
                         ihovlast = k
                         hx0 = x0
                         hx1 = x1
                         hpy = py
                      end if
                   end if
                end do
                i = j + 1
             end do
             call mo_diag_flusharrows(dl,nbuf,bufok,bufx,bufy,bufn,ic == 1,harr,coloc,&
                0.3_c_float*fontsize%x)
          end do

          ! the hovered level, marked as the table marks a hovered row.
          ! Its geometry was kept when the hit test picked it, rather than
          ! worked out again here from a second copy of the slot layout
          if (ihover > 0) &
             call mo_diag_markrect(dl,hx0,hx1,hpy,colhov,1._c_float)

          call ImDrawList_PopClipRect(dl)

          ! Selection: latched on the press and committed on the
          ! release, and only if the mouse neither moved nor was
          ! double-clicked in between. Selecting an orbital starts a grid
          ! sampling that can take seconds, so a stray selection at the
          ! start of a pan is expensive; and a double click over a level
          ! is ImPlot's fit-the-whole-spectrum gesture, whose second
          ! click would otherwise toggle off the orbital the first one
          ! had just selected and sampled
          if (inplot .and. ihover > 0 .and. igIsMouseClicked(ImGuiMouseButton_Left,.false._c_bool)) then
             w%mo_diag%ipress = ihover
             w%mo_diag%moved = .false.
             w%mo_diag%dbl = igIsMouseDoubleClicked(ImGuiMouseButton_Left)
          end if
          if (w%mo_diag%ipress /= 0) then
             if (igIsMouseDragging(ImGuiMouseButton_Left,-1._c_float)) w%mo_diag%moved = .true.
             if (igIsMouseReleased(ImGuiMouseButton_Left)) then
                if (.not.w%mo_diag%moved .and. .not.w%mo_diag%dbl .and. ihover == w%mo_diag%ipress) then
                   if (ihovlast > ihover) then
                      ! a merged pile: zoom into it rather than select
                      pad = 0.1d0 * (w%mo_diag%ene(ihovlast) - w%mo_diag%ene(ihover))
                      if (pad <= 0d0) pad = 0.05d0 * max(abs(w%mo_diag%ene(ihover)),1d0)
                      w%mo_diag%fitr(1) = w%mo_diag%ene(ihover) - pad
                      w%mo_diag%fitr(2) = w%mo_diag%ene(ihovlast) + pad
                      w%mo_diag%ifit = 3
                   elseif (w%mo_selected == w%mo_diag%imo(ihover)) then
                      w%mo_selected = 0
                      changed = .true.
                   else
                      w%mo_selected = w%mo_diag%imo(ihover)
                      w%mo_scrollto = w%mo_selected
                      changed = .true.
                   end if
                end if
                w%mo_diag%ipress = 0
                w%mo_diag%moved = .false.
                w%mo_diag%dbl = .false.
             end if
          end if
       end if
       call ipEndPlot()
    end if

    ! Which band is which, written under the plot: these name the bands
    ! the way axis tick labels name categories, and down here they
    ! cannot end up drawn over a level that happens to sit at the top of
    ! the view.
    if (nch > 1 .and. haveplot) then
       call igGetCursorScreenPos(p0)
       do ic = 1, nch
          strtxt = trim(mo_spin_name(mo_diag_ispin(nch,ic))) // c_null_char
          call igCalcTextSize(sztxt,c_loc(strtxt),c_null_ptr,.false._c_bool,-1._c_float)
          p0%x = pmin%x + (real(ic,c_float) - 0.5_c_float) * bandw - 0.5_c_float * sztxt%x
          call ImDrawList_AddText_Vec2(igGetWindowDrawList(),p0,&
             igGetColorU32_Col(ImGuiCol_Text,0.6_c_float),c_loc(strtxt),c_null_ptr)
       end do
       call igDummy(ImVec2(0._c_float,fontsize%y))
    end if

    ! what the hovered level is
    if (ihover > 0) then
       if (ihovlast > ihover) then
          strtxt = string(ihovlast-ihover+1) // " orbitals, " //&
             string(w%mo_diag%ene(ihover),'f',decimal=digits) // " to " //&
             string(w%mo_diag%ene(ihovlast),'f',decimal=digits) //&
             " (click to zoom in)" // c_null_char
       else
          imo = w%mo_diag%imo(ihover)
          call wfn%get_mo_info(imo,label,occup=occup,ener=ener)
          strtxt = "Orbital " // wfn%get_mo_name(imo)
          if (len_trim(label) > 0) strtxt = strtxt // " (" // trim(label) // ")"
          strtxt = strtxt // ", occ. " // string(occup,'f',decimal=2) //&
             ", " // string(ener*unitfactor,'f',decimal=digits) // c_null_char
       end if
       call igBeginTooltip()
       call igTextUnformatted(c_loc(strtxt),c_null_ptr)
       call igEndTooltip()
    end if

  end subroutine mo_draw_diagram

  !> Build the level list behind the energy-level diagram of the
  !> molecular-orbitals window w: one entry per orbital of the
  !> wavefunction wfn of field ifield in system isys, with its energy in
  !> the display units (unitfactor), sorted by energy within each spin
  !> channel. Returns at once if the list is still valid, so this runs
  !> once per system, per reload of the field, and per change of the
  !> energy unit -- never every frame.
  subroutine mo_diagram_build(w,wfn,isys,ifield,unitfactor)
    use systems, only: sys
    use tools, only: mergesort
    use wfn_private, only: molwfn
    type(window), intent(inout) :: w
    type(molwfn), intent(in) :: wfn
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    real*8, intent(in) :: unitfactor

    integer :: ic, ispin, nsp, i, j, imo, ntot, i0
    integer, allocatable :: iperm(:), itmp(:)
    real*8, allocatable :: etmp(:), otmp(:)
    character(len=:), allocatable :: label
    real*8 :: occup, ener

    ! nothing to do while the list still describes the field on screen
    if (w%mo_diag%isys == isys .and. w%mo_diag%gen == sys(isys)%fieldgen .and.&
       w%mo_diag%ifield == ifield .and. w%mo_diag%ieneunit == w%mo_ieneunit) return

    w%mo_diag = mo_diagram_state()
    w%mo_diag%isys = isys
    w%mo_diag%gen = sys(isys)%fieldgen
    w%mo_diag%ifield = ifield
    w%mo_diag%ieneunit = w%mo_ieneunit
    w%mo_diag%nch = wfn%get_mo_nchannels()

    ! room for every orbital of every channel
    ntot = 0
    do ic = 1, w%mo_diag%nch
       ntot = ntot + wfn%get_mo_nspin(mo_diag_ispin(w%mo_diag%nch,ic))
    end do
    if (ntot == 0) return
    allocate(w%mo_diag%imo(ntot),w%mo_diag%ene(ntot),w%mo_diag%occ(ntot))

    i0 = 0
    do ic = 1, w%mo_diag%nch
       w%mo_diag%ich0(ic) = i0 + 1
       ispin = mo_diag_ispin(w%mo_diag%nch,ic)
       nsp = wfn%get_mo_nspin(ispin)
       if (nsp == 0) cycle

       ! collect the channel
       allocate(itmp(nsp),etmp(nsp),otmp(nsp),iperm(nsp))
       j = 0
       do i = 1, nsp
          imo = wfn%get_mo_index(ispin,i)
          if (imo == 0) cycle
          call wfn%get_mo_info(imo,label,occup=occup,ener=ener)
          j = j + 1
          itmp(j) = imo
          etmp(j) = ener * unitfactor
          otmp(j) = occup
       end do

       ! sort it by energy: the file order is not guaranteed to be
       ! sorted, and the level grouping downstream needs it to be
       if (j > 0) then
          do i = 1, j
             iperm(i) = i
          end do
          call mergesort(etmp,iperm,1,j)
          do i = 1, j
             w%mo_diag%imo(i0+i) = itmp(iperm(i))
             w%mo_diag%ene(i0+i) = etmp(iperm(i))
             w%mo_diag%occ(i0+i) = otmp(iperm(i))
          end do

          ! the frontier of this channel: the highest occupied and the
          ! lowest empty level. Reading them off the occupations rather
          ! than off the orbital count covers fractional occupations too
          do i = i0+1, i0+j
             if (w%mo_diag%occ(i) > mo_diag_occtol) w%mo_diag%homo(ic) = w%mo_diag%ene(i)
          end do
          do i = i0+1, i0+j
             if (w%mo_diag%occ(i) <= mo_diag_occtol) then
                w%mo_diag%lumo(ic) = w%mo_diag%ene(i)
                exit
             end if
          end do
       end if

       i0 = i0 + j
       deallocate(itmp,etmp,otmp,iperm)
    end do
    w%mo_diag%ich0(w%mo_diag%nch+1) = i0 + 1
    w%mo_diag%n = i0
    if (i0 > 0) then
       w%mo_diag%erange(1) = minval(w%mo_diag%ene(1:i0))
       w%mo_diag%erange(2) = maxval(w%mo_diag%ene(1:i0))
    end if

  end subroutine mo_diagram_build

  !> Spin channel ic (1 or 2) of a wavefunction with nch channels, as a
  !> wfn_spin_* selector: the whole packed list when there is only one.
  function mo_diag_ispin(nch,ic) result(ispin)
    use wfn_private, only: wfn_spin_all
    integer, intent(in) :: nch
    integer, intent(in) :: ic
    integer :: ispin

    if (nch > 1) then
       ispin = ic
    else
       ispin = wfn_spin_all
    end if

  end function mo_diag_ispin

  !> Last element of the ascending array arr(lo:hi) that falls strictly
  !> below val, or lo-1 if there is none.
  function mo_diag_lower(arr,lo,hi,val) result(k)
    real(c_double), intent(in) :: arr(:)
    integer, intent(in) :: lo
    integer, intent(in) :: hi
    real*8, intent(in) :: val
    integer :: k

    integer :: ia, ib, im

    k = lo - 1
    if (hi < lo) return
    if (arr(lo) >= val) return
    if (arr(hi) < val) then
       k = hi
       return
    end if
    ia = lo
    ib = hi
    do while (ib - ia > 1)
       im = (ia + ib) / 2
       if (arr(im) < val) then
          ia = im
       else
          ib = im
       end if
    end do
    k = ia

  end function mo_diag_lower

  !> The energy range the diagram should show, answering the pending
  !> request in w%mo_diag%ifit: the region around the HOMO/LUMO gap
  !> (1), every level (2), or a range asked for by a drill-in or by the
  !> table moving the selection (3).
  subroutine mo_diagram_range(w,ylo,yhi)
    type(window), intent(in) :: w
    real(c_double), intent(out) :: ylo
    real(c_double), intent(out) :: yhi

    integer :: ic, lo, hi, k
    logical :: haveocc, havevir
    real*8 :: pad

    associate (d => w%mo_diag)
      if (d%ifit == 3) then
         ! taken as given: whoever asked for it did its own padding
         ylo = d%fitr(1)
         yhi = d%fitr(2)
         return
      end if

      if (d%ifit == 2) then
         ylo = d%erange(1)
         yhi = d%erange(2)
      else
         ! The frontier: a few levels either side of the gap. Each
         ! channel proposes its own window and the narrowest is taken,
         ! not their union -- a channel with only a handful of orbitals
         ! reaches its own core within a few levels of its HOMO, and
         ! letting it set the range would squeeze the whole valence
         ! region into a few pixels.
         haveocc = .false.
         havevir = .false.
         ylo = -huge(1d0)
         yhi = huge(1d0)
         do ic = 1, d%nch
            lo = d%ich0(ic)
            hi = d%ich0(ic+1) - 1
            if (hi < lo) cycle
            if (d%homo(ic) < huge(1d0)) then
               k = max(mo_diag_lower(d%ene,lo,hi,d%homo(ic)) + 1 - mo_diag_nshow,lo)
               ylo = max(ylo,d%ene(k))
               haveocc = .true.
            end if
            if (d%lumo(ic) < huge(1d0)) then
               k = min(mo_diag_lower(d%ene,lo,hi,d%lumo(ic)) + mo_diag_nshow,hi)
               yhi = min(yhi,d%ene(k))
               havevir = .true.
            end if
         end do
         ! with no level on one side of the gap, that side runs to the
         ! end of the spectrum: there is nothing beyond it to leave out
         if (.not.havevir) yhi = d%erange(2)
         if (.not.haveocc) ylo = d%erange(1)
         ! the windows of two channels may not overlap at all
         if (ylo > yhi) then
            ylo = d%erange(1)
            yhi = d%erange(2)
         end if
      end if

      ! a little air, and something to show for a single level
      pad = 0.05d0 * (yhi - ylo)
      if (pad <= 0d0) pad = max(0.05d0 * abs(yhi),0.05d0)
      ylo = ylo - pad
      yhi = yhi + pad
    end associate

  end subroutine mo_diagram_range

  !> Bring the orbital selected in the table into the diagram's view,
  !> keeping the zoom the user has set. Does nothing while the level is
  !> already on screen, so this cannot fight a range the user chose.
  subroutine mo_diagram_seek(w)
    type(window), intent(inout) :: w

    integer :: i
    real*8 :: span

    if (w%mo_selected <= 0 .or. w%mo_diag%n == 0) return
    do i = 1, w%mo_diag%n
       if (w%mo_diag%imo(i) /= w%mo_selected) cycle
       if (w%mo_diag%ene(i) >= w%mo_diag%ylim(1) .and. w%mo_diag%ene(i) <= w%mo_diag%ylim(2)) return
       span = w%mo_diag%ylim(2) - w%mo_diag%ylim(1)
       if (span <= 0d0) return
       w%mo_diag%fitr(1) = w%mo_diag%ene(i) - 0.5d0 * span
       w%mo_diag%fitr(2) = w%mo_diag%ene(i) + 0.5d0 * span
       w%mo_diag%ifit = 3
       return
    end do

  end subroutine mo_diagram_seek

  !> Widths of the columns of the molecular-orbital table for spin
  !> channel ispin of wavefunction wfn -- Id, cached-grid mark, Label,
  !> Spin, Occ., Energy, zero for the ones this wavefunction does not
  !> show -- and the width of the whole table.
  !>
  !> Every column is given the width of the widest entry it can ever
  !> show. Left to auto-fit, the columns whose content varies from row to
  !> row (Id and Label) resize to whatever the clipper has made visible,
  !> so they, and everything drawn to their right, shift whenever the
  !> visible set changes, which reads as the table trembling while it is
  !> scrolled. A column costs the width asked for plus the cell padding
  !> on either side, and the table also keeps a pixel per border and the
  !> width of its scroll bar: budget the total short and the last column
  !> is the one that gets clipped.
  subroutine mo_table_widths(wfn,ispin,wid,wtot)
    use gui_main, only: g
    use wfn_private, only: molwfn, wfn_uhf, wfn_rohf, wfn_spin_all
    use utils, only: iw_calcwidth
    use tools_io, only: string
    type(molwfn), intent(in) :: wfn
    integer, intent(in) :: ispin
    real(c_float), intent(out) :: wid(6)
    real(c_float), intent(out) :: wtot

    integer :: nrow, nocc, ndig, nlab, ncol
    real(c_float) :: cellpad

    nrow = wfn%get_mo_nspin(ispin,nocc)
    cellpad = 2._c_float * g%Style%CellPadding%x

    if (ispin == wfn_spin_all) then
       ndig = len_trim(string(max(wfn%nmoall,1)))
    else
       ndig = len_trim(string(max(nrow,1)))
    end if
    nlab = 5 + len_trim(string(max(max(nocc-1,nrow-nocc-1),1)))

    wid = 0._c_float
    wid(1) = iw_calcwidth(max(ndig,2),0) + cellpad
    wid(2) = 0.75_c_float * igGetFrameHeight() + cellpad
    wid(3) = iw_calcwidth(max(nlab,5),0) + cellpad
    if (ispin == wfn_spin_all .and. (wfn%wfntyp == wfn_uhf .or. wfn%wfntyp == wfn_rohf)) &
       wid(4) = iw_calcwidth(5,0) + cellpad
    wid(5) = iw_calcwidth(mo_occ_len,0) + cellpad
    if (wfn%hasene) wid(6) = iw_calcwidth(max(mo_ene_len,6),0) + cellpad

    ncol = count(wid > 0._c_float)
    wtot = sum(wid) + real(ncol,c_float) * cellpad + real(ncol+1,c_float) + g%Style%ScrollbarSize

  end subroutine mo_table_widths

  !> The button that puts the orbital table back on the HOMO/LUMO
  !> boundary. The one-shot centering that runs when the window opens is
  !> the only other thing that ever takes the table there, so once it has
  !> been scrolled away this is the way back -- and re-arming that flag is
  !> the whole implementation.
  subroutine mo_frontier_button(w,ttshown)
    use utils, only: iw_button, iw_tooltip
    type(window), intent(inout) :: w
    logical, intent(inout) :: ttshown

    if (iw_button("Frontier##mofrontier")) w%mo_scrolled = .false.
    call iw_tooltip("Scroll the table back to the HOMO/LUMO boundary",ttshown)

  end subroutine mo_frontier_button

  !> Show str on window w until it has been up long enough to read. Set
  !> bad for a failure (shown in the danger color) and clear it for a
  !> report of something that worked.
  subroutine mo_message(w,str,bad)
    use interfaces_glfw, only: glfwGetTime
    type(window), intent(inout) :: w
    character(len=*), intent(in) :: str
    logical, intent(in) :: bad

    w%errmsg = str
    w%mo_msgbad = bad
    w%timelast_errmsg = glfwGetTime()

  end subroutine mo_message

  !> Whether the cost measurement in c still describes field ifield of
  !> system isys: same system, same field, same field-set generation, and
  !> the geometry has not moved since. The system has to be part of the
  !> key -- two freshly loaded systems share a field generation and a
  !> geometry stamp, so without it a measurement made on one is taken to
  !> hold for the other, and the window silently stops coarsening.
  module function mo_cost_matches(c,isys,ifield)
    use systems, only: sys, sysc
    class(mo_cost_state), intent(in) :: c
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    logical :: mo_cost_matches

    mo_cost_matches = (c%isys == isys .and. c%ifield == ifield .and.&
       c%gen == sys(isys)%fieldgen .and. c%timegeom == sysc(isys)%timelastchange_geometry)

  end function mo_cost_matches

  !> Seconds to sample one orbital of the field window w is showing at
  !> quality ilevel (0 = automatic, else a named level), and the grid n
  !> that quality resolves to. Uses the box and the per-point cost the
  !> window has already measured, so it prices a setting without
  !> committing to it.
  function mo_quality_cost(w,isys,ilevel,n) result(secs)
    use representations, only: iso_grid_size, iso_defaultlevel, iso_npts_custom_min,&
       iso_npts_custom_max
    use param, only: bohrtoa
    type(window), intent(in) :: w
    integer, intent(in) :: isys
    integer, intent(in) :: ilevel
    integer, intent(out) :: n(3)
    real*8 :: secs

    integer :: i
    real*8 :: alen(3)

    if (ilevel == 0) then
       ! the box edge lengths, in angstrom, are the box itself
       do i = 1, 3
          alen(i) = norm2(w%mo_cost%box(:,i)) * bohrtoa
       end do
       n = iso_grid_size(isys,iso_defaultlevel,box=w%mo_cost%box,&
          ptsang=mo_auto_ptsang(w%mo_cost%secs,alen))
       ! within the bounds a custom grid obeys, since that is what the
       ! automatic quality is recorded as
       n = min(max(n,iso_npts_custom_min),iso_npts_custom_max)
    else
       n = iso_grid_size(isys,ilevel,box=w%mo_cost%box)
    end if
    secs = w%mo_cost%secs * product(real(n,8))

  end function mo_quality_cost

  !> Points per angstrom the automatic quality uses to sample one orbital
  !> over a box of edge lengths alen (angstrom), given costest, the
  !> measured cost of one sample point in seconds.
  !>
  !> The automatic quality only ever COARSENS: it is the user's default
  !> quality (the iso_defaultlevel preference), dropped to whatever the
  !> time budget allows when that level would take too long. Letting the
  !> budget also raise the resolution would make a small molecule spend
  !> the whole budget on a grid far finer than its isosurface needs --
  !> slower than the plain default, for a case that was never the
  !> problem. How coarse it can get is bounded from below by
  !> iso_grid_size's own per-axis floor, so on a pathologically slow
  !> field the budget can still be exceeded; the combo says so.
  function mo_auto_ptsang(costest,alen) result(ppa)
    use representations, only: iso_level_ptsang, iso_defaultlevel
    real*8, intent(in) :: costest
    real*8, intent(in) :: alen(3)
    real*8 :: ppa

    real*8 :: vol

    ppa = iso_level_ptsang(iso_defaultlevel)
    if (costest <= 0d0) return
    vol = alen(1) * alen(2) * alen(3)
    if (vol <= 0d0) return
    ppa = min((mo_time_budget / costest / vol)**(1d0/3d0),ppa)

  end function mo_auto_ptsang

  !> Measure the cost of one sample point of an MO of field ifield of
  !> system isys, over the box representation r samples with n points per
  !> axis, and keep it on window w. Measured once per field and geometry:
  !> the benchmark itself costs about a tenth of a second, and the answer
  !> only changes when the wavefunction or the box does.
  subroutine mo_measure_cost(w,isys,ifield,r,n)
    use systems, only: sys, sysc
    use representations, only: representation, iso_estimate_cost
    type(window), intent(inout) :: w
    integer, intent(in) :: isys
    integer, intent(in) :: ifield
    type(representation), intent(in) :: r
    integer, intent(in) :: n(3)

    ! keyed without looking at the result: a field the benchmark cannot
    ! measure must be remembered as such, or it is retried every frame
    if (w%mo_cost%matches(isys,ifield)) return

    ! benchmark what the sampling loop actually evaluates: one orbital,
    ! not the density, which costs several times more per point
    w%mo_cost%secs = iso_estimate_cost(isys,ifield,r%iso%iregion,r%iso%rgn_x,n,&
       r%iso%mo_request())
    w%mo_cost%isys = isys
    w%mo_cost%ifield = ifield
    w%mo_cost%gen = sys(isys)%fieldgen
    w%mo_cost%timegeom = sysc(isys)%timelastchange_geometry

  end subroutine mo_measure_cost

  !> The rectangle that marks a level bar running from x0 to x1 at
  !> height y as selected or hovered.
  subroutine mo_diag_markrect(dl,x0,x1,y,col,thick)
    type(c_ptr), intent(in) :: dl
    real(c_float), intent(in) :: x0
    real(c_float), intent(in) :: x1
    real(c_float), intent(in) :: y
    integer(c_int), intent(in) :: col
    real(c_float), intent(in) :: thick

    real(c_float) :: dy

    dy = 0.5_c_float * mo_diag_selthick + 2._c_float
    call ImDrawList_AddRect(dl,ImVec2(x0-2._c_float,y-dy),ImVec2(x1+2._c_float,y+dy),&
       col,2._c_float,0_c_int,thick)

  end subroutine mo_diag_markrect

  !> Draw the n electron arrows held back for one bunch of levels, or
  !> none of them if the bunch could not be laid out in full. A bunch
  !> that shows arrows on some of its levels and not on others reads as
  !> if the others were empty, which is worse than showing none: when the
  !> levels are crunched together the whole bunch goes without, and the
  !> occupations are still one hover away.
  subroutine mo_diag_flusharrows(dl,n,ok,x,y,nel,up,h,col,dx)
    type(c_ptr), intent(in) :: dl
    integer, intent(in) :: n
    logical, intent(in) :: ok
    real(c_float), intent(in) :: x(:)
    real(c_float), intent(in) :: y(:)
    integer, intent(in) :: nel(:)
    logical, intent(in) :: up
    real(c_float), intent(in) :: h
    integer(c_int), intent(in) :: col
    real(c_float), intent(in) :: dx

    integer :: i

    if (.not.ok) return
    do i = 1, n
       if (nel(i) == 2) then
          call mo_diag_arrow(dl,x(i)-dx,y(i),h,.true.,col)
          call mo_diag_arrow(dl,x(i)+dx,y(i),h,.false.,col)
       else
          call mo_diag_arrow(dl,x(i),y(i),h,up,col)
       end if
    end do

  end subroutine mo_diag_flusharrows

  !> One electron arrow of height h centered on (x,y), pointing up or
  !> down. Screen y grows downward, so an up arrow's head has the
  !> smaller y.
  subroutine mo_diag_arrow(dl,x,y,h,up,col)
    type(c_ptr), intent(in) :: dl
    real(c_float), intent(in) :: x
    real(c_float), intent(in) :: y
    real(c_float), intent(in) :: h
    logical, intent(in) :: up
    integer(c_int), intent(in) :: col

    real(c_float) :: ytail, yhead, whead, dir
    type(ImVec2) :: pa, pb, q1, q2

    if (up) then
       ytail = y + 0.5_c_float * h
       yhead = y - 0.5_c_float * h
    else
       ytail = y - 0.5_c_float * h
       yhead = y + 0.5_c_float * h
    end if
    dir = sign(1._c_float,ytail - yhead)
    whead = 0.20_c_float * h

    pa%x = x
    pa%y = ytail
    pb%x = x
    pb%y = yhead
    call ImDrawList_AddLine(dl,pa,pb,col,max(1._c_float,0.07_c_float*h))
    q1%x = x - whead
    q1%y = yhead + 2._c_float * whead * dir
    q2%x = x + whead
    q2%y = q1%y
    call ImDrawList_AddTriangleFilled(dl,pb,q1,q2,col)

  end subroutine mo_diag_arrow

end submodule mo
