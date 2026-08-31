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
  real(c_float), parameter :: mo_sep_rgb(3) = (/0.85_c_float,0.15_c_float,0.15_c_float/) ! HOMO/LUMO separator color

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
    use systems, only: sysc, sys, sys_init, ok_system, reload_field_with_virtuals
    use representations, only: representation, reptype_isosurface, repflavor_isosurface,&
       iso_isoval_mo, iso_defaultlevel, iso_level_optstr, iso_alpha_def, iso_rgb_palette,&
       iso_grid_size, iso_region_to_box
    use wfn_private, only: wfn_rhf, wfn_uhf, wfn_rohf
    use types, only: id_mo_id, field_evaluation_avail, fieldeval_category_mo
    use utils, only: iw_text, iw_button, iw_tooltip, iw_calcheight, iw_combo_simple,&
       iw_dragfloat_real8, iw_coloredit, iw_close_event, iw_setpos_bottomright,&
       iw_table_column
    use tools_io, only: string, ioj_right
    use param, only: hartoev
    class(window), intent(inout), target :: w

    logical(c_bool) :: selected
    logical :: doquit, goodsys, mo_ok, goodparent, syschanged, changed, found, ldum, hasspin
    integer :: isys, iview, iref, i, j, nrow, jsep, jscroll, itrep, irep, ihomo, spin, ncol, digits
    integer :: icid, iccache, iclabel, icspin, icocc, icene
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, strl
    character(len=:), allocatable :: label
    type(ImVec2) :: sz0, szero, sz1, p0, p1
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    real*8 :: occup, ener, unitfactor, alpha8
    type(field_evaluation_avail) :: av

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window acts on the system shown in its anchor view
    goodparent = w%anchor(iview,isys,syschanged)

    ! initialize state
    if (w%firstpass) then
       w%mo_ieneunit = 0
       w%mo_ilevel = iso_defaultlevel
       w%mo_isoval = iso_isoval_mo
       w%mo_rgb(:,1) = iso_rgb_palette(:,1) ! positive lobe: blue
       w%mo_rgb(:,2) = iso_rgb_palette(:,4) ! negative lobe: red
       w%mo_alpha = iso_alpha_def
    end if
    if (w%firstpass .or. syschanged) then
       w%errmsg = "" ! the error, if any, was about the previous system
       w%mo_selected = 0
       w%mo_scrolled = .false.
       w%mo_cache = mo_cache_state() ! the grids belong to the previous system
    end if

    ! initialize
    szero%x = 0
    szero%y = 0
    changed = .false.
    itrep = 0
    doquit = .not.goodparent
    if (.not.doquit) doquit = .not.associated(win(iview)%sc)

    ! can the reference field of this system provide molecular orbitals?
    goodsys = ok_system(isys,sys_init)
    mo_ok = goodsys
    iref = 0
    if (mo_ok) then
       iref = sys(isys)%iref
       mo_ok = sys(isys)%goodfield(iref)
    end if
    if (mo_ok) then
       call sys(isys)%f(iref)%eval_avail(av)
       mo_ok = av%avail(fieldeval_category_mo)
    end if
    if (mo_ok) mo_ok = associated(win(iview)%sc)

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
       if (.not.mo_ok) &
          call iw_text("The reference field cannot provide molecular orbitals",danger=.true.,wrap=.true.)
    end if

    ! maybe the error message
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

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
             ! the occupied orbitals are unchanged: keep their grids
             if (len_trim(w%errmsg) == 0) call mo_cache_keep_after_reload(isys,iref)
          end if
          call iw_tooltip("Re-read the wavefunction file for the reference field,&
             & this time including the virtual (unoccupied) orbitals",ttshown)
       elseif (.not.sys(isys)%f(iref)%wfn%hasvirtual) then
          call iw_text("No virtuals available")
       end if

       associate (wfn => sys(isys)%f(iref)%wfn)
         ! drop a stale selection
         if (w%mo_selected > wfn%nmoall) w%mo_selected = 0
         ihomo = wfn%nmoocc

         ! orbital table header and energy units
         call iw_text("Orbitals",highlight=.true.,alignframe=.true.)
         if (wfn%hasene) then
            call iw_combo_simple("##moeneunit","Hartree" // c_null_char // "eV" // c_null_char,&
               w%mo_ieneunit,sameline=.true.)
            call iw_tooltip("Units for the orbital energies",ttshown)
         end if
         if (w%mo_ieneunit == 0) then
            unitfactor = 1d0
            digits = 4
         else
            unitfactor = hartoev
            digits = 3
         end if

         ! the orbitals table; the spin and energy columns appear only
         ! when the wavefunction provides them
         hasspin = (wfn%wfntyp == wfn_uhf .or. wfn%wfntyp == wfn_rohf)
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

         ! The table runs in orbital order, lowest first. When there are
         ! virtuals, a separator row marks the HOMO/LUMO gap and the
         ! one-shot centering puts that boundary in the middle of the
         ! table; without virtuals there is no gap, and the HOMO (the
         ! last row) is centered instead.
         nrow = wfn%nmoall
         jsep = 0
         if (wfn%nmoall > wfn%nmoocc) then
            jsep = wfn%nmoocc + 1
            nrow = nrow + 1
            jscroll = jsep
         else
            jscroll = ihomo
         end if

         flags = ImGuiTableFlags_None
         flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
         flags = ior(flags,ImGuiTableFlags_RowBg)
         flags = ior(flags,ImGuiTableFlags_Borders)
         flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
         flags = ior(flags,ImGuiTableFlags_ScrollY)
         str1 = "##tablemolecularorbitals" // c_null_char
         sz0%x = 0._c_float
         sz0%y = iw_calcheight(10,0,.false.)
         if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
            ! header setup
            call iw_table_column("Id",id=icid,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("",id=iccache,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("Label",id=iclabel,flags=ImGuiTableColumnFlags_WidthFixed)
            if (hasspin) &
               call iw_table_column("Spin",id=icspin,flags=ImGuiTableColumnFlags_WidthFixed)
            call iw_table_column("Occ.",id=icocc,flags=ImGuiTableColumnFlags_WidthFixed)
            if (wfn%hasene) &
               call iw_table_column("Energy",id=icene,flags=ImGuiTableColumnFlags_WidthFixed)
            call igTableSetupScrollFreeze(0, 1) ! top row always visible

            ! draw the header
            call igTableHeadersRow()
            call igTableSetColumnWidthAutoAll(igGetCurrentTable())

            ! draw the rows through a clipper (the wavefunction may
            ! have thousands of MOs); while the one-shot centering on
            ! the HOMO is pending, force its row in so it can anchor
            ! the scroll
            clipper = ImGuiListClipper_ImGuiListClipper()
            call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)
            if (.not.w%mo_scrolled) &
               call ImGuiListClipper_ForceDisplayRangeByIndices(clipper,jscroll-1,jscroll)
            do while (ImGuiListClipper_Step(clipper))
               call c_f_pointer(clipper,clipper_f)
               do j = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                  ! the HOMO/LUMO gap: a thick line across the table
                  if (j == jsep) then
                     call igTableNextRow(ImGuiTableRowFlags_None,mo_sep_frac*igGetFrameHeight())
                     call igTableSetBgColor(ImGuiTableBgTarget_RowBg0,&
                        igGetColorU32_Vec4(ImVec4(mo_sep_rgb(1),mo_sep_rgb(2),mo_sep_rgb(3),&
                        1._c_float)),-1_c_int)
                     if (.not.w%mo_scrolled .and. j == jscroll) then
                        if (igTableSetColumnIndex(icid)) call igSetScrollHereY(0.5_c_float)
                        w%mo_scrolled = .true.
                     end if
                     cycle
                  end if

                  ! display row to orbital (the rows below the separator
                  ! are shifted by it)
                  i = j
                  if (jsep > 0 .and. j > jsep) i = i - 1

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
                     if (.not.w%mo_scrolled .and. j == jscroll) then
                        call igSetScrollHereY(0.5_c_float)
                        w%mo_scrolled = .true.
                     end if

                     ! text
                     call iw_text(string(i),sameline=.true.)
                  end if

                  ! whether this orbital's sampling grid is in the cache:
                  ! a dot as tall as it is wide, filling the height of the
                  ! row. Drawn on top of the row selectable (emitted in the
                  ! id column, before this one) so it survives the
                  ! selection highlight; the dummy reserves the width the
                  ! column is auto-sized to
                  if (igTableSetColumnIndex(iccache)) then
                     sz1%x = igGetFrameHeight()
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
                     if (igTableSetColumnIndex(icspin)) then
                        if (spin == 1) then
                           call iw_text("alpha")
                        elseif (spin == 2) then
                           call iw_text("beta")
                        else
                           call iw_text("both")
                        end if
                     end if
                  end if

                  ! occupation
                  if (igTableSetColumnIndex(icocc)) &
                     call iw_text(string(occup,'f',length=6,decimal=2,justify=ioj_right))

                  ! energy
                  if (wfn%hasene) then
                     if (igTableSetColumnIndex(icene)) &
                        call iw_text(string(ener*unitfactor,'f',length=12,decimal=digits,justify=ioj_right))
                  end if
               end do ! clipper range
            end do ! clipper step
            call ImGuiListClipper_End(clipper)
            call ImGuiListClipper_destroy(clipper)
            call igEndTable()
         end if ! igBeginTable

         ! display options
         call iw_text("Display",highlight=.true.)

         ! grid quality
         call iw_text("Quality",alignframe=.true.)
         call iw_combo_simple("##molevel",iso_level_optstr,w%mo_ilevel,startsatone=.true.,&
            sameline=.true.,changed=ldum)
         changed = changed .or. ldum
         call iw_tooltip("Density of the grid used to calculate the MO isosurface",ttshown)

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
                    if (r%iso%ilevel /= w%mo_ilevel .or. .not.r%iso%isgenerated(isys)) then
                       r%iso%ilevel = w%mo_ilevel
                       call commit_grid(win(iview)%sc%reptrans(itrep))
                    end if
                    win(iview)%sc%forcebuildlists = .true.
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
               if (len_trim(label) > 0) then
                  win(iview)%sc%rep(irep)%name = "MO " // string(w%mo_selected) // " (" // trim(label) // ")"
               else
                  win(iview)%sc%rep(irep)%name = "MO " // string(w%mo_selected)
               end if
               win(iview)%sc%forcebuildlists = .true.
            end if
         end if
         call iw_tooltip("Copy the displayed orbital into a permanent isosurface object,&
            & editable from the Objects menu of the view",ttshown)
       end associate
    end if ! mo_ok

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
    subroutine commit_grid(r)
      type(representation), intent(inout) :: r

      integer :: n(3)
      real*8 :: box(3,0:3)
      logical :: okbox

      call iso_region_to_box(isys,r%iso%iregion,r%iso%rgn_x,box,okbox)
      if (.not.okbox) return
      n = iso_grid_size(isys,r%iso%ilevel,r%iso%nptscustom,box=box)
      call r%iso%apply_grid(n,r%iso%iregion,r%iso%rgn_x)

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

end submodule mo
