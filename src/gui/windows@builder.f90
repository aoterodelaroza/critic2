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
    use systems, only: sysc, ok_system, sys_init
    use utils, only: iw_text, iw_button, iw_tooltip
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_PICKATOM_SELECT,&
       BIND_PICKATOM_ALT, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use interfaces_glfw, only: glfwGetTime
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: isys, iview, icel, imode
    logical :: doquit, goodparent, havesys, ok

    logical, save :: ttshown = .false. ! tooltip flag

    real*8, parameter :: stale_gap = 1d0 ! discard clicks older than this (s)

    ! do we have a good parent window (a view)? A hidden view is still
    ! good: the builder survives the hide and the mode can be released
    iview = w%idparent
    goodparent = iview > 0 .and. iview <= nwin
    if (goodparent) goodparent = win(iview)%isinit
    if (goodparent) goodparent = (win(iview)%type == wintype_view)
    doquit = .not.goodparent

    ! initialize the state on first pass
    if (w%firstpass) then
       w%builder_active = .false.
       w%builder_isys = 0
       w%builder_time = 0d0
    end if

    ! the builder operates on the system shown in the parent view
    isys = 0
    havesys = .false.
    if (goodparent) then
       isys = win(iview)%isys
       havesys = ok_system(isys,sys_init)
    end if

    ! handle an active change-valence mode commanded to the parent view
    if (w%builder_active) then
       ok = goodparent
       if (ok) ok = win(iview)%viewmode == vm_builder .and.&
          win(iview)%vmdata%owner == w%id .and.&
          win(iview)%isys == w%builder_isys
       if (.not.ok) then
          ! the view is gone, the mode was exited (any key, mode combo),
          ! another window took over the pick, or the view shows another
          ! system: release the forced mode
          call builder_stop()
       elseif (win(iview)%vmdata%idx(1) > 0) then
          ! an atom click was delivered: apply the edit (main pick = add a
          ! hydrogen, alternate pick = remove one) and stay in the mode.
          ! Discard stale clicks instead: those delivered while the builder
          ! was not polling (hidden dock tab, collapsed window) or with a
          ! geometry change after the last click-free poll, because the
          ! stored cell index may then refer to a different atom.
          icel = win(iview)%vmdata%idx(1)
          imode = win(iview)%vmdata%flag
          win(iview)%vmdata%idx = 0
          if (glfwGetTime() - w%builder_time < stale_gap .and.&
             sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) &
             call sysc(w%builder_isys)%change_valence(icel,imode)
       else
          ! no click pending: stamp the reference time for the stale-click guard
          w%builder_time = glfwGetTime()
       end if
    end if

    ! header: the system this builder operates on
    if (havesys) then
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)
    end if

    ! the valence section
    call iw_text("Valence",highlight=.true.)
    if (iw_button("Change",disabled=.not.havesys)) &
       call builder_toggle()
    call iw_tooltip("Add ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//") or remove ("//&
       trim(get_bind_keyname(BIND_PICKATOM_ALT))//") hydrogens on the clicked atoms in the"//&
       " view, repositioning the terminal substituents (click again to stop)",ttshown)

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
    ! Stop the change-valence mode: release the parent view (if alive)
    ! and clear the mode state.
    subroutine builder_stop()

      if (goodparent) call win(iview)%viewmode_release_forced(w%id,vm_builder)
      w%builder_active = .false.
      w%builder_isys = 0
    end subroutine builder_stop

    ! Toggle the change-valence mode: start it on the parent view, or
    ! stop it if it is already active.
    subroutine builder_toggle()

      if (w%builder_active) then
         call builder_stop()
      else
         w%builder_active = .true.
         w%builder_isys = isys
         w%builder_time = glfwGetTime()
         call win(iview)%viewmode_set_forced(vm_builder,&
            "Add ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//&
            ") or remove ("//trim(get_bind_keyname(BIND_PICKATOM_ALT))//&
            ") hydrogens (press any key to exit)...",w%id)
      end if
    end subroutine builder_toggle
  end subroutine draw_builder

end submodule builder
