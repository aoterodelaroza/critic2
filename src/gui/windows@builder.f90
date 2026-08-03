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
    use gui_main, only: ColorHighlightEditDistScene
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple, iw_dragfloat_real8
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_PICKATOM_SELECT,&
       BIND_PICKATOM_ALT, BIND_RECALC_BONDS, BIND_NAV_MEASURE, BIND_EDITDISTANCE,&
       BIND_REOPEN, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use interfaces_glfw, only: glfwGetTime
    use tools_io, only: string
    use param, only: bohrtoa
    class(window), intent(inout), target :: w

    integer :: isys, iview, icel, imode, nsel, iside, ibold, ibnew
    logical :: doquit, goodparent, havesys, ok, lcommit, ldum
    real*8 :: dist
    character(len=:), allocatable :: stropt

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

    ! initialize the state on first pass (editdist_pending is NOT reset:
    ! the view sets it right before this window's first draw)
    if (w%firstpass) then
       w%errmsg = ""
       w%builder_vm = 0
       w%builder_isys = 0
       w%builder_time = 0d0
       w%editdist_active = .false.
       w%editdist_isys = 0
       w%editdist_idx = 0
       w%editdist_dirty = .false.
       w%editdist_time = 0d0
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
       ! drop any edit-distance session without rebuilding (the geometry is
       ! discarded) and treat the system as not ready for the rest of this
       ! frame: the re-read deallocates sys(isys)%c until the init thread runs
       w%editdist_dirty = .false.
       if (w%editdist_active) call editdist_stop()
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

    ! edit distance section
    if (w%editdist_pending) then
       w%editdist_pending = .false.
       if (.not.w%editdist_active) call editdist_start()
    end if
    if (w%editdist_active) then
       ok = havesys
       if (ok) ok = (isys == w%editdist_isys)
       if (ok) ok = all(w%editdist_idx(1,:) >= 1) .and. all(w%editdist_idx(1,:) <= sys(isys)%c%ncel)
       if (ok) ok = .not.sysc(isys)%md_run
       ! any external geometry or bond edit invalidates the latched atoms
       ! and masked fragments (our own edits re-stamp editdist_time)
       if (ok) ok = sysc(isys)%timelastchange_geometry <= w%editdist_time
       if (ok) ok = sysc(isys)%timelastchange_rebond <= w%editdist_time
       if (.not.ok) then
          call editdist_stop()
       else
          w%editdist_time = glfwGetTime()
       end if
    end if

    ! number of measure-selected atoms in the parent view
    nsel = 0
    if (goodparent) then
       if (associated(win(iview)%sc)) nsel = win(iview)%sc%nmsel
    end if

    call iw_text("Edit distance",highlight=.true.)
    if (.not.w%editdist_active) &
       call iw_text("Select two atoms in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//&
       "), then press Edit distance or "//trim(get_bind_keyname(BIND_EDITDISTANCE)))
    ok = w%editdist_active .or. (havesys .and. nsel == 2)
    if (iw_button("Edit distance",disabled=.not.ok)) then
       w%errmsg = ""
       if (w%editdist_active) then
          call editdist_stop()
       else
          call editdist_start()
       end if
    end if
    call iw_tooltip("Edit the distance and bond between the two selected atoms ("//&
       trim(get_bind_keyname(BIND_EDITDISTANCE))//"); press again to stop",ttshown)

    if (w%editdist_active) then
       ! the two atoms
       call iw_text("(1) "//trim(sys(isys)%c%at(sys(isys)%c%atcel(w%editdist_idx(1,1))%idx)%name)//&
          string(w%editdist_idx(1,1)))
       call iw_text("(2) "//trim(sys(isys)%c%at(sys(isys)%c%atcel(w%editdist_idx(1,2))%idx)%name)//&
          string(w%editdist_idx(1,2)),sameline=.true.)

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
          if (w%editdist_fragok(iside)) stropt = stropt // "Translate fragment"//c_null_char
          call iw_combo_simple("Atom "//string(iside)//"##editdistmove"//string(iside),&
             stropt,w%editdist_imove(iside))
          call iw_tooltip("What moves when the distance changes: nothing (fixed), the atom,"//&
             " or the whole fragment attached to it",ttshown)
       end do

       ! distance drag-float (bohr internally, shown in angstrom, no upper bound)
       dist = norm2(editdist_pos(2) - editdist_pos(1))
       if (all(w%editdist_imove == 0)) call igBeginDisabled(.true._c_bool)
       ldum = iw_dragfloat_real8("Distance (Å)##editdistdrag",x1=dist,speed=0.01d0,&
          min=0.1d0,scale=bohrtoa,decimal=3,notlive=.true.,committed=lcommit,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       if (all(w%editdist_imove == 0)) then
          call igEndDisabled()
       else
          if (ldum) call editdist_apply(dist)
          if (lcommit) call editdist_commit()
       end if
       call iw_tooltip("Distance between the two atoms (Å)",ttshown)

       ! apply button
       if (iw_button("Apply",danger=.true.)) then
          w%errmsg = ""
          call editdist_stop()
       end if
       call iw_tooltip("Keep the current distance and release the two atoms",ttshown)
    end if

    ! transient highlight of the two latched atoms
    if (w%editdist_active) &
       call sysc(w%editdist_isys)%highlight_atoms(.true.,w%editdist_idx(1,:),&
       atlisttype_ncel_frac,spread(ColorHighlightEditDistScene,2,2))

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
    ! Start the edit-distance session from the parent view's two
    ! measure-selected atoms: latch the atoms and the system, compute
    ! the masked fragments of both sides, and set the defaults.
    subroutine editdist_start()
      integer :: ia, ib
      logical :: disc
      integer, allocatable :: ifr1(:), ifr2(:)

      if (.not.havesys) return
      if (.not.associated(win(iview)%sc)) return
      if (win(iview)%sc%nmsel /= 2) return
      ia = win(iview)%sc%msel(1,1)
      ib = win(iview)%sc%msel(1,2)
      if (ia < 1 .or. ia > sys(isys)%c%ncel .or. ib < 1 .or. ib > sys(isys)%c%ncel) return
      if (ia == ib) then
         w%errmsg = "The selected atoms are periodic images of the same atom"
         return
      end if
      if (sysc(isys)%md_run) then
         w%errmsg = "Edit distance is not available while dynamics is running"
         return
      end if

      ! latch the two atoms (cell atom + lattice vector) and the system,
      ! and clear the measurement that provided them
      w%editdist_idx(:,1) = win(iview)%sc%msel(1:4,1)
      w%editdist_idx(:,2) = win(iview)%sc%msel(1:4,2)
      w%editdist_isys = isys
      win(iview)%sc%nmsel = 0
      win(iview)%sc%msel = 0

      ! masked fragments: connected components ignoring all bonds
      ! between the two atoms. The fragment option is unavailable if
      ! the component is periodic or contains the other atom.
      ! masked_fragment seeds the list with the atom itself, so the
      ! translate-atom mode is the one-atom head of the same list.
      call sys(isys)%c%masked_fragment(ia,ia,ib,w%editdist_nfrag(1),ifr1,disc)
      w%editdist_fragok(1) = disc .and. .not.any(ifr1(1:w%editdist_nfrag(1)) == ib)
      call sys(isys)%c%masked_fragment(ib,ia,ib,w%editdist_nfrag(2),ifr2,disc)
      w%editdist_fragok(2) = disc .and. .not.any(ifr2(1:w%editdist_nfrag(2)) == ia)
      if (allocated(w%editdist_frag)) deallocate(w%editdist_frag)
      allocate(w%editdist_frag(maxval(w%editdist_nfrag),2))
      w%editdist_frag = 0
      w%editdist_frag(1:w%editdist_nfrag(1),1) = ifr1(1:w%editdist_nfrag(1))
      w%editdist_frag(1:w%editdist_nfrag(2),2) = ifr2(1:w%editdist_nfrag(2))

      ! defaults: atom 1 fixed; atom 2 translates its fragment if
      ! available, else the atom
      w%editdist_imove(1) = 0
      w%editdist_imove(2) = merge(2,1,w%editdist_fragok(2))

      w%editdist_dirty = .false.
      w%editdist_time = glfwGetTime()
      w%editdist_active = .true.
    end subroutine editdist_start

    ! Stop the edit-distance session, rebuilding the structure first if
    ! in-place moves were applied but not committed (editdist_start
    ! re-seeds the rest of the session state).
    subroutine editdist_stop()

      if (w%editdist_dirty) then
         if (ok_system(w%editdist_isys,sys_init)) then
            call sys(w%editdist_isys)%c%rebuild_after_move(copybonding=.true.)
            if (goodparent) then
               if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
            end if
            call sysc(w%editdist_isys)%post_event(lastchange_geometry)
         end if
         w%editdist_dirty = .false.
      end if
      w%editdist_active = .false.
      w%editdist_isys = 0
      w%editdist_idx = 0
    end subroutine editdist_stop

    ! Cartesian position (bohr) of the latched image of atom iside.
    function editdist_pos(iside) result(r)
      integer, intent(in) :: iside
      real*8 :: r(3)

      r = sys(isys)%c%atcel(w%editdist_idx(1,iside))%r +&
         sys(isys)%c%x2c(dble(w%editdist_idx(2:4,iside)))
    end function editdist_pos

    ! Index in the bond combo (0 = none, else bondorder(idx)) of the
    ! current bond between the two latched atoms.
    function editdist_bondidx() result(idx)
      integer :: idx
      integer :: ia, j, k, lv(3)

      idx = 0
      ia = w%editdist_idx(1,1)
      lv = w%editdist_idx(2:4,2) - w%editdist_idx(2:4,1)
      if (.not.allocated(sys(isys)%c%nstar)) return
      do j = 1, sys(isys)%c%nstar(ia)%ncon
         if (sys(isys)%c%nstar(ia)%idcon(j) == w%editdist_idx(1,2) .and.&
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

      lv = w%editdist_idx(2:4,2) - w%editdist_idx(2:4,1)
      if (inew == 0) then
         call sysc(isys)%remove_bond(w%editdist_idx(1,1),w%editdist_idx(1,2),lv)
      elseif (iold == 0) then
         call sysc(isys)%add_bond(w%editdist_idx(1,1),w%editdist_idx(1,2),lv,bondorder(inew))
      else
         call sysc(isys)%set_bond_order(w%editdist_idx(1,1),w%editdist_idx(1,2),lv,bondorder(inew))
      end if
      w%editdist_time = glfwGetTime()
    end subroutine editdist_setbond

    ! Move the atoms so the distance between the two latched images
    ! becomes dnew (bohr), displacing each side along the bond axis
    ! according to its move mode (in-place update; the rebuild happens
    ! at editdist_commit).
    subroutine editdist_apply(dnew)
      real*8, intent(in) :: dnew

      integer :: i, k, n
      real*8 :: ra(3), rb(3), u(3), dcur, delta, dshift(2)
      real*8, allocatable :: rnew(:,:)

      ra = editdist_pos(1)
      rb = editdist_pos(2)
      dcur = norm2(rb - ra)
      if (dcur < eps_dzero) return
      u = (rb - ra) / dcur
      delta = dnew - dcur

      ! split the displacement among the movable sides (atom 1 moves
      ! against the bond axis, atom 2 along it)
      if (w%editdist_imove(1) > 0 .and. w%editdist_imove(2) > 0) then
         dshift = (/-delta/2d0, delta/2d0/)
      elseif (w%editdist_imove(1) > 0) then
         dshift = (/-delta, 0d0/)
      elseif (w%editdist_imove(2) > 0) then
         dshift = (/0d0, delta/)
      else
         return
      end if

      ! build the new positions and apply them in place; the moving set
      ! is the masked-fragment list (translate fragment) or its one-atom
      ! head (translate atom)
      allocate(rnew(3,sys(isys)%c%ncel))
      do i = 1, sys(isys)%c%ncel
         rnew(:,i) = sys(isys)%c%atcel(i)%r
      end do
      do i = 1, 2
         if (w%editdist_imove(i) == 0) cycle
         n = merge(w%editdist_nfrag(i),1,w%editdist_imove(i) == 2)
         do k = 1, n
            rnew(:,w%editdist_frag(k,i)) = rnew(:,w%editdist_frag(k,i)) + dshift(i) * u
         end do
      end do
      call sys(isys)%c%update_positions(rnew)
      w%editdist_dirty = .true.
      if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
      call sysc(isys)%post_event(lastchange_geometry)
      ! stamp after the post so our own edit does not end the session
      w%editdist_time = glfwGetTime()
    end subroutine editdist_apply

    ! Commit the edit-distance moves: rebuild the structure from the
    ! moved positions (keeping the current bonds) and re-latch the
    ! lattice vectors (crystal atoms may be rewrapped by the rebuild).
    subroutine editdist_commit()
      integer :: i
      real*8 :: rab(3,2)

      do i = 1, 2
         rab(:,i) = editdist_pos(i)
      end do
      call sys(isys)%c%rebuild_after_move(copybonding=.true.)
      if (.not.sys(isys)%c%ismolecule) then
         do i = 1, 2
            w%editdist_idx(2:4,i) = nint(sys(isys)%c%c2x(rab(:,i)) -&
               sys(isys)%c%atcel(w%editdist_idx(1,i))%x)
         end do
      end if
      w%editdist_dirty = .false.
      if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
      call sysc(isys)%post_event(lastchange_geometry)
      w%editdist_time = glfwGetTime()
    end subroutine editdist_commit

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

end submodule builder
