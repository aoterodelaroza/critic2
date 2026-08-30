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
contains

  !> Draw the molecular orbitals window. The window lists the MOs of
  !> the reference field of the system in its anchor view; selecting
  !> one shows it as a transient +/- isosurface pair in that view
  !> (kept alive by re-arming it every frame, and reaped when the
  !> selection or the window goes away). The Create Isosurface Object
  !> button copies the transient into a permanent representation.
  module subroutine draw_mo(w)
    use systems, only: sysc, sys, sys_init, ok_system
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
    integer :: isys, iview, iref, i, itrep, irep, ihomo, spin, ncol, digits
    integer :: icid, iclabel, icspin, icocc, icene
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, strl
    character(len=:), allocatable :: label
    type(ImVec2) :: sz0, szero
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
       w%mo_selected = 0
       w%mo_scrolled = .false.
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
         if (wfn%hasvirtual) then
            call iw_text(", " // string(wfn%nmoall-wfn%nmoocc) // " virtual",sameline_nospace=.true.)
         else
            call iw_text(" (no virtuals available)",sameline_nospace=.true.)
         end if

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
         iclabel = 1
         ncol = 2
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
            call ImGuiListClipper_Begin(clipper,wfn%nmoall,-1._c_float)
            if (.not.w%mo_scrolled) &
               call ImGuiListClipper_ForceDisplayRangeByIndices(clipper,ihomo-1,ihomo)
            do while (ImGuiListClipper_Step(clipper))
               call c_f_pointer(clipper,clipper_f)
               do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
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

                     ! center the table on the HOMO the first time it is drawn
                     if (.not.w%mo_scrolled .and. i == ihomo) then
                        call igSetScrollHereY(0.5_c_float)
                        w%mo_scrolled = .true.
                     end if

                     ! text
                     call iw_text(string(i),sameline=.true.)
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

end submodule mo
