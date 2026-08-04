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

! Routines for the builder window.
submodule (windows) builder
  use interfaces_cimgui
  implicit none

contains

  !> Draw the builder window. The builder is tied to its parent view
  !> (w%idparent) and operates on the system shown in it.
  module subroutine draw_builder(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_geometry,&
       reread_system_from_file, atlisttype_ncel_frac
    use dynamics, only: md_relax
    use energy, only: ff_backend_applicable, ff_backend_default
    use gui_main, only: ColorHighlightEditDistScene
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple, iw_dragfloat_real8
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_PICKATOM_SELECT,&
       BIND_PICKATOM_ALT, BIND_RECALC_BONDS, BIND_NAV_MEASURE, BIND_EDIT_D_A_PHI,&
       BIND_REOPEN, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG, BIND_CANCEL
    use interfaces_glfw, only: glfwGetTime
    use tools_io, only: string
    use tools_math, only: cross, axisangle2mat
    use param, only: bohrtoa, pi
    class(window), intent(inout), target :: w

    integer :: isys, iview, icel, imode, nsel, iside, ibold, ibnew
    logical :: doquit, goodparent, havesys, ok, lcommit, ldum, relaxing
    real*8 :: dist, ang
    character(len=:), allocatable :: stropt, errmsg

    logical, save :: ttshown = .false. ! tooltip flag

    real*8, parameter :: stale_gap = 1d0 ! discard clicks older than this (s)
    real*8, parameter :: eps_dzero = 1d-10 ! degenerate-distance threshold (bohr)
    ! bond order for each bond-combo index >= 1 (single, double, triple, dashed, aromatic)
    integer, parameter :: bondorder(5) = (/1,2,3,0,-1/)

    ! do we have a good parent window (a view)? A hidden view is still
    ! good: the builder survives the hide and the mode can be released
    iview = w%idparent
    goodparent = iview > 0 .and. iview <= nwin
    if (goodparent) goodparent = win(iview)%isinit
    if (goodparent) goodparent = (win(iview)%type == wintype_view)
    doquit = .not.goodparent

    ! initialize the state on first pass (edit_pending is NOT reset:
    ! the view sets it right before this window's first draw)
    if (w%firstpass) then
       w%errmsg = ""
       w%builder_vm = 0
       w%builder_isys = 0
       w%builder_time = 0d0
       w%edit_kind = 0
       w%edit_isys = 0
       w%edit_idx = 0
       w%edit_dirty = .false.
       w%edit_time = 0d0
    end if

    ! the builder operates on the system shown in the parent view
    isys = 0
    havesys = .false.
    if (goodparent) then
       isys = win(iview)%isys
       havesys = ok_system(isys,sys_init)
    end if

    ! handle an active builder mode commanded to the parent view
    if (w%builder_vm /= 0) then
       ok = goodparent
       if (ok) ok = win(iview)%viewmode == w%builder_vm .and.&
          win(iview)%vmdata%owner == w%id .and.&
          win(iview)%isys == w%builder_isys
       if (.not.ok) then
          ! The view is gone, the mode was exited (any key, mode combo),
          ! another window took over the pick, or the view shows another
          ! system: release the forced mode
          call builder_stop()
       elseif (win(iview)%vmdata%idx(1) > 0) then
          ! An atom click was delivered: apply the edit and stay in
          ! the mode. Discard stale clicks: those delivered while the
          ! builder was not polling.
          icel = win(iview)%vmdata%idx(1)
          imode = win(iview)%vmdata%flag
          win(iview)%vmdata%idx = 0
          if (glfwGetTime() - w%builder_time < stale_gap .and.&
             sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) then
             if (w%builder_vm == vm_builder_valence) then
                ! change valence: main pick = add a hydrogen, alternate
                ! pick = remove one
                call sysc(w%builder_isys)%change_valence(icel,imode)
             elseif (w%builder_vm == vm_builder_remove .and. imode == 1) then
                ! remove atoms (main pick only): the atom and its
                ! terminal hydrogens
                call sysc(w%builder_isys)%remove_atom_hydrogens(icel)
             end if
             ! hold the camera of the view the user is clicking in through
             ! the rebuild (the edit routines post the geometry event)
             if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
          end if
       else
          ! no click pending: stamp the reference time for the stale-click guard
          w%builder_time = glfwGetTime()
       end if
    end if

    ! header: the system this builder operates on
    call iw_text("System",highlight=.true.)
    if (havesys) &
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)
    if (iw_button("Restore",danger=.true.,disabled=.not.havesys)) then
       w%errmsg = ""
       ! drop any edit session without rebuilding (the geometry is
       ! discarded) and treat the system as not ready for the rest of this
       ! frame: the re-read deallocates sys(isys)%c until the init thread runs
       w%edit_dirty = .false.
       call w%edit_stop()
       call reread_system_from_file(isys)
       havesys = .false.
    end if
    call iw_tooltip("Restore the system to the original geometry it had when it was"//&
       " first opened ("//trim(get_bind_keyname(BIND_REOPEN))//")",ttshown)

    ! the atoms section
    call iw_text("Atoms",highlight=.true.)
    if (iw_button("Remove atoms",disabled=.not.havesys)) then
       w%errmsg = ""
       call builder_toggle(vm_builder_remove)
    end if
    call iw_tooltip("Remove ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//") the clicked"//&
       " atoms in the view, along with their terminal hydrogens",ttshown)

    ! the valence section
    call iw_text("Valence",highlight=.true.)
    if (iw_button("Change",disabled=.not.havesys)) then
       w%errmsg = ""
       call builder_toggle(vm_builder_valence)
    end if
    call iw_tooltip("Add ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//") or remove ("//&
       trim(get_bind_keyname(BIND_PICKATOM_ALT))//") hydrogens on the clicked atoms in the"//&
       " view, repositioning the terminal substituents",ttshown)

    ! the bonds section
    call iw_text("Bonds",highlight=.true.)
    if (iw_button("Rebond",disabled=.not.havesys)) then
       w%errmsg = ""
       call sysc(isys)%rebond()
    end if
    call iw_tooltip("Recompute the bond connectivity for this system ("//&
       trim(get_bind_keyname(BIND_RECALC_BONDS))//")",ttshown)

    ! number of measure-selected atoms in the parent view
    nsel = 0
    if (goodparent) then
       if (associated(win(iview)%sc)) nsel = win(iview)%sc%nmsel
    end if

    ! the keybinding request toggles an edit session: two selected
    ! atoms toggle edit distance, three edit angle, four edit dihedral,
    ! and otherwise deactivate the active session, if any
    if (w%edit_pending .or. (w%focused() .and. is_bind_event(BIND_EDIT_D_A_PHI))) then
       w%edit_pending = .false.
       if (nsel >= 2 .and. nsel <= 4) then
          call edit_toggle(nsel)
       elseif (w%edit_kind /= 0) then
          call w%edit_stop()
       end if
    end if

    ! per-frame guards on the active edit session
    if (w%edit_kind /= 0) then
       if (edit_session_ok()) then
          w%edit_time = glfwGetTime()
       else
          call w%edit_stop()
       end if
    end if

    ! edit distance section
    call iw_text("Edit distance",highlight=.true.)
    if (w%edit_kind /= 2) &
       call iw_text("Select two atoms in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//&
       "), then press Edit distance or "//trim(get_bind_keyname(BIND_EDIT_D_A_PHI)))
    ok = (w%edit_kind == 2) .or. (havesys .and. nsel == 2)
    if (iw_button("Edit distance",disabled=.not.ok)) then
       w%errmsg = ""
       call edit_toggle(2)
    end if
    call iw_tooltip("Edit the distance and bond between the two selected atoms ("//&
       trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//"); press again to stop",ttshown)

    if (w%edit_kind == 2) then
       ! the two atoms
       call edit_atom_labels()

       ! bond type combo
       ibold = editdist_bondidx()
       ibnew = ibold
       call iw_combo_simple("Bond##editdistbond","None"//c_null_char//"Single"//c_null_char//&
          "Double"//c_null_char//"Triple"//c_null_char//"Dashed"//c_null_char//&
          "Aromatic"//c_null_char,ibnew)
       call iw_tooltip("Type of the bond between the two atoms",ttshown)
       if (ibnew /= ibold) call editdist_setbond(ibold,ibnew)

       ! per-atom move mode combos
       do iside = 1, 2
          stropt = "Fixed"//c_null_char//"Translate atom"//c_null_char
          if (w%edit_fragok(iside)) stropt = stropt // "Translate fragment"//c_null_char
          call iw_combo_simple("Atom "//string(iside)//"##editdistmove"//string(iside),&
             stropt,w%edit_imove(iside))
          call iw_tooltip("What moves when the distance changes: nothing (fixed), the atom,"//&
             " or the whole fragment attached to it",ttshown)
       end do

       ! distance drag-float (bohr internally, shown in angstrom, no upper bound)
       dist = norm2(edit_pos(2) - edit_pos(1))
       if (all(w%edit_imove == 0)) call igBeginDisabled(.true._c_bool)
       ldum = iw_dragfloat_real8("Distance (Å)##editdistdrag",x1=dist,speed=0.01d0,&
          min=0.1d0,scale=bohrtoa,decimal=3,notlive=.true.,committed=lcommit,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       if (all(w%edit_imove == 0)) then
          call igEndDisabled()
       else
          if (ldum) call editdist_apply(dist)
          if (lcommit) call edit_commit()
       end if
       call iw_tooltip("Distance between the two atoms (Å)",ttshown)

       ! apply button
       if (iw_button("Apply",danger=.true.)) then
          w%errmsg = ""
          call w%edit_stop()
       end if
       call iw_tooltip("Keep the current distance and release the two atoms",ttshown)
    end if

    ! edit angle section
    call iw_text("Edit angle",highlight=.true.)
    if (w%edit_kind /= 3) &
       call iw_text("Select three atoms in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//&
       "), then press Edit angle or "//trim(get_bind_keyname(BIND_EDIT_D_A_PHI)))
    ok = (w%edit_kind == 3) .or. (havesys .and. nsel == 3)
    if (iw_button("Edit angle",disabled=.not.ok)) then
       w%errmsg = ""
       call edit_toggle(3)
    end if
    call iw_tooltip("Edit the angle between the three selected atoms, vertex at the"//&
       " second atom ("//trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//"); press again to stop",ttshown)

    if (w%edit_kind == 3) then
       ! the three atoms
       call edit_atom_labels()

       ! central atom combo; while the central atom moves, the
       ! terminals are held fixed (rotation needs a fixed vertex)
       stropt = "Fixed"//c_null_char//"Translate atom"//c_null_char
       if (w%edit_fragok(2)) stropt = stropt // "Translate group"//c_null_char
       call iw_combo_simple("Atom 2 (center)##editangmove2",stropt,w%edit_imove(2))
       call iw_tooltip("What moves when the angle changes: nothing (fixed), the central"//&
          " atom, or the whole group attached to it",ttshown)
       if (w%edit_imove(2) > 0) then
          w%edit_imove(1) = 0
          w%edit_imove(3) = 0
       end if

       ! terminal atom combos (rotation available only with a fixed center)
       do iside = 1, 3, 2
          stropt = "Fixed"//c_null_char//"Rotate atom"//c_null_char
          if (w%edit_fragok(iside)) stropt = stropt // "Rotate group"//c_null_char
          if (w%edit_imove(2) > 0) call igBeginDisabled(.true._c_bool)
          call iw_combo_simple("Atom "//string(iside)//"##editangmove"//string(iside),&
             stropt,w%edit_imove(iside))
          if (w%edit_imove(2) > 0) call igEndDisabled()
          call iw_tooltip("What moves when the angle changes: nothing (fixed), the"//&
             " terminal atom, or the whole group attached to it, rotating about the"//&
             " central atom",ttshown)
       end do

       ! angle drag-float (degrees)
       ang = editang_angat(edit_pos(2)) * 180d0 / pi
       if (all(w%edit_imove == 0)) call igBeginDisabled(.true._c_bool)
       ldum = iw_dragfloat_real8("Angle (°)##editangdrag",x1=ang,speed=0.2d0,&
          min=0d0,max=180d0,decimal=2,notlive=.true.,committed=lcommit,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       if (all(w%edit_imove == 0)) then
          call igEndDisabled()
       else
          if (ldum) call editang_apply(ang)
          if (lcommit) call edit_commit()
       end if
       call iw_tooltip("Angle between the three atoms (degrees)",ttshown)

       ! apply button
       if (iw_button("Apply##editang",danger=.true.)) then
          w%errmsg = ""
          call w%edit_stop()
       end if
       call iw_tooltip("Keep the current angle and release the three atoms",ttshown)
    end if

    ! edit dihedral section
    call iw_text("Edit dihedral",highlight=.true.)
    if (w%edit_kind /= 4) &
       call iw_text("Select four atoms in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//&
       "), then press Edit dihedral or "//trim(get_bind_keyname(BIND_EDIT_D_A_PHI)))
    ok = (w%edit_kind == 4) .or. (havesys .and. nsel == 4)
    if (iw_button("Edit dihedral",disabled=.not.ok)) then
       w%errmsg = ""
       call edit_toggle(4)
    end if
    call iw_tooltip("Edit the dihedral angle of the four selected atoms about the 2-3"//&
       " bond ("//trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//"); press again to stop",ttshown)

    if (w%edit_kind == 4) then
       ! the four atoms
       call edit_atom_labels()

       ! terminal atom combos (atoms 2 and 3 always stay fixed)
       do iside = 1, 4, 3
          stropt = "Fixed"//c_null_char//"Rotate atom"//c_null_char
          if (w%edit_fragok(iside)) stropt = stropt // "Rotate group (fix 2-3)"//c_null_char
          ! fragok(2)/(3) = validity of the half adjacent to terminal 1/4
          if (w%edit_fragok(merge(2,3,iside == 1))) &
             stropt = stropt // "Rotate group (move 2-3)"//c_null_char
          call iw_combo_simple("Atom "//string(iside)//"##editdihmove"//string(iside),&
             stropt,w%edit_imove(iside))
          call iw_tooltip("What rotates about the 2-3 axis when the dihedral changes:"//&
             " nothing (fixed), the terminal atom, the group attached to it (the"//&
             " substituents of atoms 2 and 3 stay fixed), or the whole half severed"//&
             " at the 2-3 bond (the environments of atoms 2 and 3 stay rigid)",ttshown)
       end do

       ! dihedral drag-float (degrees)
       ang = editdih_val() * 180d0 / pi
       if (all(w%edit_imove == 0)) call igBeginDisabled(.true._c_bool)
       ldum = iw_dragfloat_real8("Dihedral (°)##editdihdrag",x1=ang,speed=0.2d0,&
          min=-180d0,max=180d0,decimal=2,notlive=.true.,committed=lcommit,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       if (all(w%edit_imove == 0)) then
          call igEndDisabled()
       else
          if (ldum) call editdih_apply(ang)
          if (lcommit) call edit_commit()
       end if
       call iw_tooltip("Dihedral angle of the four atoms (degrees)",ttshown)

       ! apply button
       if (iw_button("Apply##editdih",danger=.true.)) then
          w%errmsg = ""
          call w%edit_stop()
       end if
       call iw_tooltip("Keep the current dihedral and release the four atoms",ttshown)
    end if

    ! transient highlight of the latched atoms
    if (w%edit_kind /= 0) &
       call sysc(w%edit_isys)%highlight_atoms(.true.,w%edit_idx(1,1:w%edit_kind),&
       atlisttype_ncel_frac,spread(ColorHighlightEditDistScene,2,w%edit_kind))

    ! the relax section
    call iw_text("Relax",highlight=.true.)
    relaxing = .false.
    if (havesys) then
       relaxing = sysc(isys)%md_run .and. sysc(isys)%md%mode == md_relax
       ! resolve the method before the button (the combo below would
       ! only snap an inapplicable backend after this frame's click)
       if (sysc(isys)%md_backend < 0 .or.&
          .not.ff_backend_applicable(sysc(isys)%md_backend,sys(isys)%c)) &
          sysc(isys)%md_backend = ff_backend_default(sys(isys)%c)
    end if
    if (.not.relaxing) then
       ok = havesys
       if (ok) ok = .not.sysc(isys)%md_run
       if (iw_button("Relax",danger=.true.,disabled=.not.ok)) then
          ! release any edit session (a running relaxation ends it anyway)
          call w%edit_stop()
          call sysc(isys)%md_start(md_relax,errmsg)
          w%errmsg = errmsg
          if (len_trim(w%errmsg) == 0) then
             if (associated(win(iview)%sc)) win(iview)%forcerender = .true.
          end if
       end if
       call iw_tooltip("Relax the geometry to the nearest energy minimum, stopping when"//&
          " the maximum force falls below the threshold",ttshown)
    else
       if (iw_button("Stop",danger=.true.)) then
          call sysc(isys)%md_stop()
          if (associated(win(iview)%sc)) win(iview)%forcerender = .true.
       end if
       call iw_tooltip("Stop the geometry relaxation ("//&
          trim(get_bind_keyname(BIND_CANCEL))//")",ttshown)
    end if

    if (havesys) then
       ! method combo, next to the button
       call draw_ff_backend_combo(isys,"##builderrelaxmethod",15)
       call iw_tooltip("Method for the calculation of energies and forces",ttshown)

       ! force convergence threshold
       ldum = iw_dragfloat_real8("Force threshold (eV/Å)##relaxfconv",x1=sysc(isys)%md%fconv,&
          speed=0.0005d0,min=0.0001d0,max=1d0,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Stop the relaxation when the maximum force falls below this value",ttshown)

       ! live status: energy and maximum force
       if (sysc(isys)%md%ready .and. sysc(isys)%md%nat > 0 .and.&
          sysc(isys)%md%mode == md_relax) then
          call iw_text("Energy: "//string(sysc(isys)%md%epot,'f',decimal=6)//" Ha")
          call iw_text("Max force: "//string(sysc(isys)%md%maxforce(),'f',decimal=4)//" eV/Å")
          if (.not.sysc(isys)%md_run .and. sysc(isys)%md%converged()) &
             call iw_text("(converged)",highlight=.true.,sameline=.true.)
       end if

       ! surface force-evaluation errors from the run
       if (allocated(sysc(isys)%md%errmsg)) then
          if (len_trim(sysc(isys)%md%errmsg) > 0) &
             call iw_text(trim(sysc(isys)%md%errmsg),danger=.true.)
       end if
    end if

    ! the symmetry section
    call iw_text("Symmetry",highlight=.true.)
    if (iw_button("Refine symmetry",disabled=.not.havesys)) then
       w%errmsg = ""
       call sysc(isys)%refine_symmetry(w%errmsg)
       ! hold the camera on the parent view's scene (the main scene when
       ! the parent is the main view, its private scene otherwise)
       if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
    end if
    call iw_tooltip("Refine the geometry having atoms occupy their symmetry positions exactly",ttshown)

    ! error message, if any
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! close button and binds
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit the window (window_end releases the forced mode if still active)
    if (doquit) &
       call w%end()

  contains
    ! Per-frame validity of the active edit session: the parent view
    ! must show the same system, the latched atoms must exist, dynamics
    ! must not be running, and no external geometry or bond edit may
    ! have occurred (our own edits re-stamp the heartbeat).
    function edit_session_ok() result(ok)
      logical :: ok

      ok = havesys
      if (ok) ok = (isys == w%edit_isys)
      if (ok) ok = all(w%edit_idx(1,1:w%edit_kind) >= 1) .and.&
         all(w%edit_idx(1,1:w%edit_kind) <= sys(isys)%c%ncel)
      if (ok) ok = .not.sysc(isys)%md_run
      if (ok) ok = sysc(isys)%timelastchange_geometry <= w%edit_time
      if (ok) ok = sysc(isys)%timelastchange_rebond <= w%edit_time
    end function edit_session_ok

    ! Toggle an edit session of the given kind (2 = distance, 3 =
    ! angle, 4 = dihedral): stop it if it is the active session, and
    ! start it from the measure-selected atoms otherwise.
    subroutine edit_toggle(ikind)
      integer, intent(in) :: ikind

      if (w%edit_kind == ikind) then
         call w%edit_stop()
      else
         call edit_start(ikind)
      end if
    end subroutine edit_toggle

    ! Start an edit session of the given kind (2 = distance, 3 =
    ! angle, 4 = dihedral) from the parent view's measure-selected
    ! atoms: latch the atoms and the system, compute the masked
    ! groups, and set the defaults.
    subroutine edit_start(ikind)
      integer, intent(in) :: ikind

      integer :: idx1(4), i, j
      character(len=:), allocatable :: kname

      if (.not.havesys) return
      if (.not.associated(win(iview)%sc)) return
      if (win(iview)%sc%nmsel /= ikind) return
      idx1 = 0
      idx1(1:ikind) = win(iview)%sc%msel(1,1:ikind)
      if (minval(idx1(1:ikind)) < 1 .or. maxval(idx1(1:ikind)) > sys(isys)%c%ncel) return
      do i = 1, ikind-1
         do j = i+1, ikind
            if (idx1(i) /= idx1(j)) cycle
            ! for dihedrals, atoms 2 and 3 never move, so two periodic
            ! images of the same atom are allowed there (chain crystals)
            if (ikind == 4 .and. i == 2 .and. j == 3) cycle
            w%errmsg = "Two of the selected atoms are periodic images of the same atom"
            return
         end do
      end do
      if (sysc(isys)%md_run) then
         if (ikind == 2) then
            kname = "distance"
         elseif (ikind == 3) then
            kname = "angle"
         else
            kname = "dihedral"
         end if
         w%errmsg = "Edit " // kname // " is not available while dynamics is running"
         return
      end if
      call edit_latch(ikind)
      call edit_latch_groups(ikind)

      ! defaults: all atoms fixed except the last, which moves its
      ! group if available, else the atom; for dihedrals, prefer the
      ! group mode that also rotates the 2-3 substituents (fragok(3)
      ! is the validity of the half adjacent to terminal 4)
      w%edit_imove = 0
      if (ikind == 4 .and. w%edit_fragok(3)) then
         w%edit_imove(4) = 3
      else
         w%edit_imove(ikind) = merge(2,1,w%edit_fragok(ikind))
      end if

      w%edit_dirty = .false.
      w%edit_time = glfwGetTime()
      w%edit_kind = ikind
    end subroutine edit_start

    ! Show the labels of the latched atoms in one row.
    subroutine edit_atom_labels()
      integer :: is

      do is = 1, w%edit_kind
         call iw_text("("//string(is)//") "//&
            trim(sys(isys)%c%at(sys(isys)%c%atcel(w%edit_idx(1,is))%idx)%name)//&
            string(w%edit_idx(1,is)),sameline=(is > 1))
      end do
    end subroutine edit_atom_labels

    ! Rotate the moving terminal atoms (or their groups) about the
    ! axis through center: the first latched atom rotates by
    ! -dang/nrot and the last by +dang/nrot, which opens the angle or
    ! dihedral by dang. In the dihedral "move 2-3" mode (imove = 3)
    ! the move set is the corresponding half, stored in the interior
    ! columns.
    subroutine edit_rotate_terminals(rnew,center,axis,dang,nrot)
      real*8, intent(inout) :: rnew(:,:)
      real*8, intent(in) :: center(3), axis(3), dang
      integer, intent(in) :: nrot

      integer :: i, k, n, icol
      real*8 :: rot(3,3)

      do i = 1, w%edit_kind, w%edit_kind - 1
         if (w%edit_imove(i) == 0) cycle
         rot = axisangle2mat(axis,merge(-1d0,1d0,i == 1) * dang / nrot)
         if (w%edit_imove(i) == 3) then
            icol = merge(2,3,i == 1)
            n = w%edit_nfrag(icol)
         else
            icol = i
            n = merge(w%edit_nfrag(i),1,w%edit_imove(i) == 2)
         end if
         do k = 1, n
            rnew(:,w%edit_frag(k,icol)) = center + matmul(rot,rnew(:,w%edit_frag(k,icol)) - center)
         end do
      end do
    end subroutine edit_rotate_terminals

    ! Latch the ikind measure-selected atoms of the parent view into
    ! edit_idx/edit_isys and clear the measurement. Stops any active
    ! session first; if its uncommitted moves rebuild the crystal, the
    ! latched lattice vectors are re-derived against the rewrapped
    ! atoms.
    subroutine edit_latch(ikind)
      integer, intent(in) :: ikind

      integer :: i, ml(4,4)
      real*8 :: rabs(3,4)

      ml = 0
      ml(:,1:ikind) = win(iview)%sc%msel(1:4,1:ikind)
      do i = 1, ikind
         rabs(:,i) = sys(isys)%c%atcel(ml(1,i))%r + sys(isys)%c%x2c(dble(ml(2:4,i)))
      end do
      call w%edit_stop()
      if (.not.sys(isys)%c%ismolecule) then
         do i = 1, ikind
            ml(2:4,i) = nint(sys(isys)%c%c2x(rabs(:,i)) - sys(isys)%c%atcel(ml(1,i))%x)
         end do
      end if
      w%edit_idx = ml
      w%edit_isys = isys
      win(iview)%sc%nmsel = 0
      win(iview)%sc%msel = 0
    end subroutine edit_latch

    ! Compute the masked groups of the ikind latched atoms: each atom
    ! severs its bonds to its neighbors in the 1-2-...-ikind chain, and
    ! its group move option is available only if the component is
    ! discrete and contains none of the other latched atoms.
    ! masked_fragment seeds each list with the atom itself, so the atom
    ! move mode is the one-atom head of the same list. For dihedrals,
    ! columns 2-3 instead hold the two halves from severing the 2-3
    ! bond: they are the move sets of the terminals' "move 2-3" modes.
    subroutine edit_latch_groups(ikind)
      integer, intent(in) :: ikind

      integer :: is, js, it, n
      logical :: okf, disc
      integer, allocatable :: ifr(:), work(:,:)

      allocate(work(sys(isys)%c%ncel,ikind))
      w%edit_nfrag = 0
      w%edit_fragok = .false.
      do is = 1, ikind
         if (ikind == 4 .and. (is == 2 .or. is == 3)) then
            ! dihedral interior columns: the two halves obtained by
            ! severing the 2-3 bond, seeded at atoms 2 and 3. Each must
            ! contain its own terminal (1 for the 2-half, 4 for the
            ! 3-half) and no atom of the other half.
            it = merge(1,4,is == 2)
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,2),&
               w%edit_idx(1,3),n,ifr,disc)
         elseif (is > 1 .and. is < ikind) then
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is-1),n,ifr,disc,w%edit_idx(1,is),w%edit_idx(1,is+1))
         elseif (is > 1) then
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is-1),n,ifr,disc)
         else
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is+1),n,ifr,disc)
         end if
         ! valid if discrete, contains the seed's terminal (trivially
         ! true when it = is: the fragment contains its own seed), and
         ! contains no other latched atom
         okf = disc .and. any(ifr(1:n) == w%edit_idx(1,it))
         do js = 1, ikind
            if (js /= is .and. js /= it) okf = okf .and. .not.any(ifr(1:n) == w%edit_idx(1,js))
         end do
         w%edit_nfrag(is) = n
         w%edit_fragok(is) = okf
         work(1:n,is) = ifr(1:n)
      end do
      if (ikind == 4) then
         ! offer the "move 2-3" mode only when the plain group mode is
         ! also available: the combo option lists are built contiguously,
         ! so entry 3 (move) must not exist without entry 2 (fix). This
         ! can suppress a valid half when the terminal's own group is
         ! invalid (e.g. a ring through the terminal-axis bond).
         if (.not.w%edit_fragok(1)) w%edit_fragok(2) = .false.
         if (.not.w%edit_fragok(4)) w%edit_fragok(3) = .false.
      end if
      if (allocated(w%edit_frag)) deallocate(w%edit_frag)
      allocate(w%edit_frag(maxval(w%edit_nfrag(1:ikind)),ikind))
      w%edit_frag = 0
      do is = 1, ikind
         w%edit_frag(1:w%edit_nfrag(is),is) = work(1:w%edit_nfrag(is),is)
      end do
    end subroutine edit_latch_groups

    ! Cartesian position (bohr) of the latched image of atom iside.
    function edit_pos(iside) result(r)
      integer, intent(in) :: iside
      real*8 :: r(3)

      r = sys(isys)%c%atcel(w%edit_idx(1,iside))%r +&
         sys(isys)%c%x2c(dble(w%edit_idx(2:4,iside)))
    end function edit_pos

    ! Index in the bond combo (0 = none, else bondorder(idx)) of the
    ! current bond between the two latched atoms.
    function editdist_bondidx() result(idx)
      integer :: idx
      integer :: ia, j, k, lv(3)

      idx = 0
      ia = w%edit_idx(1,1)
      lv = w%edit_idx(2:4,2) - w%edit_idx(2:4,1)
      if (.not.allocated(sys(isys)%c%nstar)) return
      do j = 1, sys(isys)%c%nstar(ia)%ncon
         if (sys(isys)%c%nstar(ia)%idcon(j) == w%edit_idx(1,2) .and.&
            all(sys(isys)%c%nstar(ia)%lcon(:,j) == lv)) then
            do k = 1, size(bondorder,1)
               if (sys(isys)%c%nstar(ia)%ordcon(j) == bondorder(k)) then
                  idx = k
                  return
               end if
            end do
            return
         end if
      end do
    end function editdist_bondidx

    ! Change the bond between the latched atoms from combo index iold
    ! to combo index inew (0 = no bond). Posts a rebond event only, so
    ! the selection survives; re-stamp the session heartbeat so our own
    ! rebond does not end the session.
    subroutine editdist_setbond(iold,inew)
      integer, intent(in) :: iold, inew

      integer :: lv(3)

      lv = w%edit_idx(2:4,2) - w%edit_idx(2:4,1)
      if (inew == 0) then
         call sysc(isys)%remove_bond(w%edit_idx(1,1),w%edit_idx(1,2),lv)
      elseif (iold == 0) then
         call sysc(isys)%add_bond(w%edit_idx(1,1),w%edit_idx(1,2),lv,bondorder(inew))
      else
         call sysc(isys)%set_bond_order(w%edit_idx(1,1),w%edit_idx(1,2),lv,bondorder(inew))
      end if
      w%edit_time = glfwGetTime()
    end subroutine editdist_setbond

    ! Move the atoms so the distance between the two latched images
    ! becomes dnew (bohr), displacing each side along the bond axis
    ! according to its move mode (in-place update; the rebuild happens
    ! at edit_commit).
    subroutine editdist_apply(dnew)
      real*8, intent(in) :: dnew

      integer :: i, k, n
      real*8 :: ra(3), rb(3), u(3), dcur, delta, dshift(2)
      real*8, allocatable :: rnew(:,:)

      ra = edit_pos(1)
      rb = edit_pos(2)
      dcur = norm2(rb - ra)
      if (dcur < eps_dzero) return
      u = (rb - ra) / dcur
      delta = dnew - dcur

      ! split the displacement among the movable sides (atom 1 moves
      ! against the bond axis, atom 2 along it)
      if (w%edit_imove(1) > 0 .and. w%edit_imove(2) > 0) then
         dshift = (/-delta/2d0, delta/2d0/)
      elseif (w%edit_imove(1) > 0) then
         dshift = (/-delta, 0d0/)
      elseif (w%edit_imove(2) > 0) then
         dshift = (/0d0, delta/)
      else
         return
      end if

      ! build the new positions and apply them in place; the moving set
      ! is the masked-fragment list (translate fragment) or its one-atom
      ! head (translate atom)
      call edit_snapshot(rnew)
      do i = 1, 2
         if (w%edit_imove(i) == 0) cycle
         n = merge(w%edit_nfrag(i),1,w%edit_imove(i) == 2)
         do k = 1, n
            rnew(:,w%edit_frag(k,i)) = rnew(:,w%edit_frag(k,i)) + dshift(i) * u
         end do
      end do
      call edit_push(rnew)
    end subroutine editdist_apply

    ! Snapshot the current atomic positions (Cartesian, bohr) into rnew.
    subroutine edit_snapshot(rnew)
      real*8, allocatable, intent(inout) :: rnew(:,:)

      integer :: i

      allocate(rnew(3,sys(isys)%c%ncel))
      do i = 1, sys(isys)%c%ncel
         rnew(:,i) = sys(isys)%c%atcel(i)%r
      end do
    end subroutine edit_snapshot

    ! Apply the moved positions in place, mark the session dirty, and
    ! notify (the rebuild happens at edit_commit).
    subroutine edit_push(rnew)
      real*8, intent(in) :: rnew(:,:)

      call sys(isys)%c%update_positions(rnew)
      w%edit_dirty = .true.
      call edit_notify()
    end subroutine edit_push

    ! Hold the camera of the parent view through the next list rebuild,
    ! post the geometry event, and re-stamp the session heartbeat
    ! (after the post, so our own edit does not end the session).
    subroutine edit_notify()

      if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
      call sysc(isys)%post_event(lastchange_geometry)
      w%edit_time = glfwGetTime()
    end subroutine edit_notify

    ! Angle (radians) at vertex position rr subtended by the two
    ! latched terminal atoms.
    function editang_angat(rr) result(a)
      real*8, intent(in) :: rr(3)
      real*8 :: a

      real*8 :: va(3), vc(3)

      va = edit_pos(1) - rr
      vc = edit_pos(3) - rr
      a = atan2(norm2(cross(va,vc)),dot_product(va,vc))
    end function editang_angat

    ! Unit vector perpendicular to v (assumed nonzero).
    function editang_perp(v) result(u)
      real*8, intent(in) :: v(3)
      real*8 :: u(3)

      if (abs(v(1)) < 0.9d0 * norm2(v)) then
         u = cross(v,(/1d0,0d0,0d0/))
      else
         u = cross(v,(/0d0,1d0,0d0/))
      end if
      u = u / norm2(u)
    end function editang_perp

    ! Find the translation (shift) of the central atom along the
    ! internal bisector that makes the angle equal to atgt (radians),
    ! with the terminal atoms fixed. Returns .false. if no such shift
    ! can be found.
    function editang_center_shift(atgt,shift) result(ok)
      real*8, intent(in) :: atgt
      real*8, intent(out) :: shift(3)
      logical :: ok

      integer :: i
      real*8 :: ra(3), rb(3), rc(3), va(3), vc(3), w0(3), pc(3), acur
      real*8 :: t0, t1, tm, f0, f1, fm

      ok = .false.
      shift = 0d0
      ra = edit_pos(1)
      rb = edit_pos(2)
      rc = edit_pos(3)
      va = ra - rb
      vc = rc - rb
      if (norm2(va) < eps_dzero .or. norm2(vc) < eps_dzero) return
      acur = editang_angat(rb)
      f0 = acur - atgt

      if (f0 < 0d0) then
         ! increase the angle: the internal bisector from the vertex
         ! crosses the 1-3 chord (angle-bisector theorem), where the
         ! angle tends to 180. Bracket between the vertex position and
         ! a point just short of the chord.
         pc = ra + (norm2(va) / (norm2(va) + norm2(vc))) * (rc - ra)
         t1 = norm2(pc - rb)
         if (t1 < eps_dzero) return
         w0 = (pc - rb) / t1
         t1 = t1 * (1d0 - 1d-10)
         f1 = editang_angat(rb + t1 * w0) - atgt
         if (f0 * f1 > 0d0) then
            ! the target is numerically unreachable (atgt = 180): land
            ! on the chord point, where the angle is closest to it
            shift = t1 * w0
            ok = .true.
            return
         end if
      else
         ! decrease the angle: move away from the chord along the
         ! bisector (any perpendicular to the arms if collinear); the
         ! angle goes to zero, so a doubling search brackets any
         ! reachable target
         w0 = va / norm2(va) + vc / norm2(vc)
         if (norm2(w0) < eps_dzero) then
            w0 = editang_perp(va)
         else
            w0 = w0 / norm2(w0)
         end if
         t1 = -0.05d0
         do i = 1, 60
            f1 = editang_angat(rb + t1 * w0) - atgt
            if (f0 * f1 <= 0d0) exit
            t1 = 2d0 * t1
         end do
         if (f0 * f1 > 0d0) return
      end if

      ! bisection between t0 = 0 and t1
      ok = .true.
      t0 = 0d0
      do i = 1, 80
         tm = 0.5d0 * (t0 + t1)
         fm = editang_angat(rb + tm * w0) - atgt
         if (f0 * fm > 0d0) then
            t0 = tm
            f0 = fm
         else
            t1 = tm
         end if
      end do
      shift = 0.5d0 * (t0 + t1) * w0
    end function editang_center_shift

    ! Move the atoms so the angle at the central atom becomes anew
    ! (degrees): translate the center along the bisector, or rotate
    ! the terminals about the center, according to the move modes
    ! (in-place update; the rebuild happens at edit_commit).
    subroutine editang_apply(anew)
      real*8, intent(in) :: anew

      integer :: k, n, nrot
      real*8 :: rb(3), va(3), vc(3), axis(3), acur, atgt, dang
      real*8 :: shift(3)
      real*8, allocatable :: rnew(:,:)

      rb = edit_pos(2)
      va = edit_pos(1) - rb
      vc = edit_pos(3) - rb
      if (norm2(va) < eps_dzero .or. norm2(vc) < eps_dzero) return
      acur = editang_angat(rb)
      atgt = anew * pi / 180d0
      if (abs(atgt - acur) < 1d-12) return
      if (w%edit_imove(2) > 0) then
         ! translate the central atom (or its group) along the internal
         ! bisector until the angle becomes atgt
         if (.not.editang_center_shift(atgt,shift)) return
      else
         nrot = count((/w%edit_imove(1),w%edit_imove(3)/) > 0)
         if (nrot == 0) return
         dang = atgt - acur
         ! rotation axis: perpendicular to the plane of the three atoms;
         ! any direction perpendicular to the arms if they are collinear
         axis = cross(va,vc)
         if (norm2(axis) < eps_dzero) axis = editang_perp(va)
      end if

      ! build the new positions and apply them in place; the moving set
      ! is the masked-group list or its one-atom head
      call edit_snapshot(rnew)
      if (w%edit_imove(2) > 0) then
         n = merge(w%edit_nfrag(2),1,w%edit_imove(2) == 2)
         do k = 1, n
            rnew(:,w%edit_frag(k,2)) = rnew(:,w%edit_frag(k,2)) + shift
         end do
      else
         call edit_rotate_terminals(rnew,rb,axis,dang,nrot)
      end if

      call edit_push(rnew)
    end subroutine editang_apply

    ! Signed dihedral angle (radians) of the four latched atoms about
    ! the 2-3 axis, with the same convention as the measurements
    ! (0 for degenerate, collinear input).
    function editdih_val() result(a)
      real*8 :: a

      a = sys(isys)%c%dihedral(&
         sys(isys)%c%c2x(edit_pos(1)),sys(isys)%c%c2x(edit_pos(2)),&
         sys(isys)%c%c2x(edit_pos(3)),sys(isys)%c%c2x(edit_pos(4)))
    end function editdih_val

    ! Move the atoms so the dihedral angle becomes anew (degrees):
    ! rotate the terminal atoms (or their groups) about the 2-3 axis
    ! according to the move modes (in-place update; the rebuild happens
    ! at edit_commit).
    subroutine editdih_apply(anew)
      real*8, intent(in) :: anew

      integer :: nrot
      real*8 :: r2(3), u(3), acur, dang, xn
      real*8, allocatable :: rnew(:,:)

      nrot = count((/w%edit_imove(1),w%edit_imove(4)/) > 0)
      if (nrot == 0) return

      ! the 2-3 axis; the dihedral is undefined if either arm is
      ! collinear with it
      r2 = edit_pos(2)
      u = edit_pos(3) - r2
      xn = norm2(u)
      if (xn < eps_dzero) return
      u = u / xn
      if (norm2(cross(edit_pos(1) - r2,u)) < eps_dzero) return
      if (norm2(cross(edit_pos(4) - edit_pos(3),u)) < eps_dzero) return

      acur = editdih_val()
      dang = anew * pi / 180d0 - acur
      if (abs(dang) < 1d-12) return

      ! build the new positions and apply them in place
      call edit_snapshot(rnew)
      call edit_rotate_terminals(rnew,r2,u,dang,nrot)
      call edit_push(rnew)
    end subroutine editdih_apply

    ! Commit the edit-session moves: rebuild the structure from the
    ! moved positions (keeping the current bonds) and re-latch the
    ! lattice vectors (crystal atoms may be rewrapped by the rebuild).
    subroutine edit_commit()
      integer :: i
      real*8 :: rl(3,4)

      do i = 1, w%edit_kind
         rl(:,i) = edit_pos(i)
      end do
      call sys(isys)%c%rebuild_after_move(copybonding=.true.)
      if (.not.sys(isys)%c%ismolecule) then
         do i = 1, w%edit_kind
            w%edit_idx(2:4,i) = nint(sys(isys)%c%c2x(rl(:,i)) -&
               sys(isys)%c%atcel(w%edit_idx(1,i))%x)
         end do
      end if
      w%edit_dirty = .false.
      call edit_notify()
    end subroutine edit_commit

    ! Stop the active builder mode: release the parent view (if alive)
    ! and clear the mode state.
    subroutine builder_stop()

      if (goodparent) call win(iview)%viewmode_release_forced(w%id,w%builder_vm)
      w%builder_vm = 0
      w%builder_isys = 0
    end subroutine builder_stop

    ! Toggle builder mode jvm (a vm_* constant): start it on the parent
    ! view, or stop it if it is already the active mode. Starting a mode
    ! while another one is active switches to the new mode.
    subroutine builder_toggle(jvm)
      integer, intent(in) :: jvm

      character(len=:), allocatable :: msg

      if (w%builder_vm == jvm) then
         call builder_stop()
      else
         if (w%builder_vm /= 0) call builder_stop()
         w%builder_vm = jvm
         w%builder_isys = isys
         w%builder_time = glfwGetTime()
         if (jvm == vm_builder_valence) then
            msg = "Add ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//&
               ") or remove ("//trim(get_bind_keyname(BIND_PICKATOM_ALT))//&
               ") hydrogens (any key exits)..."
         else
            msg = "Remove ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//&
               ") the atoms and their hydrogens"//&
               " (any key exits)..."
         end if
         call win(iview)%viewmode_set_forced(jvm,msg,w%id)
      end if
    end subroutine builder_toggle
  end subroutine draw_builder

  !> Stop the active edit session (distance, angle, or dihedral) of a builder
  !> window: rebuild the structure first if in-place moves were applied
  !> but not committed, and release the latched atoms (the start
  !> routines re-seed the rest of the session state).
  module subroutine edit_stop(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_geometry
    class(window), intent(inout) :: w

    integer :: iview

    if (w%edit_dirty) then
       if (ok_system(w%edit_isys,sys_init)) then
          call sys(w%edit_isys)%c%rebuild_after_move(copybonding=.true.)
          iview = w%idparent
          if (iview >= 1 .and. iview <= nwin) then
             if (win(iview)%isinit) then
                if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
             end if
          end if
          call sysc(w%edit_isys)%post_event(lastchange_geometry)
       end if
       w%edit_dirty = .false.
    end if
    w%edit_kind = 0
    w%edit_isys = 0
    w%edit_idx = 0

  end subroutine edit_stop

end submodule builder
