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

  !> Draw the builder window
  module subroutine draw_builder(w)
    use systems, only: sysc, ok_system, sys_init
    use utils, only: iw_text, iw_button, iw_tooltip
    use interfaces_glfw, only: glfwGetTime
    class(window), intent(inout), target :: w

    integer :: isys, iview, ivm, icel
    logical :: havesys, okview, ok

    logical, save :: ttshown = .false. ! tooltip flag

    real*8, parameter :: stale_gap = 1d0 ! discard clicks older than this (s)

    ! initialize the state on first pass
    if (w%firstpass) then
       w%builder_mode = 0
       w%builder_iview = 0
       w%builder_isys = 0
       w%builder_time = 0d0
    end if

    ! the builder operates on the system shown in the main view
    isys = win(iwin_view)%view_selected
    havesys = win(iwin_view)%isinit .and. win(iwin_view)%isopen
    if (havesys) havesys = ok_system(isys,sys_init)

    ! handle an active valence mode commanded to a view window
    if (w%builder_mode > 0) then
       ivm = vm_builder_inc
       if (w%builder_mode == 2) ivm = vm_builder_dec
       iview = w%builder_iview
       okview = (iview >= 1 .and. iview <= nwin)
       if (okview) okview = win(iview)%isinit .and. win(iview)%isopen .and.&
          win(iview)%type == wintype_view
       ok = okview
       if (ok) ok = win(iview)%viewmode == ivm .and.&
          win(iview)%vmdata%owner == w%id .and.&
          win(iview)%view_selected == w%builder_isys
       if (.not.ok) then
          ! the view is gone, the mode was exited (any key, mode combo),
          ! another window took over the pick, or the view shows another
          ! system: release the forced mode
          if (okview) call win(iview)%viewmode_release_forced(w%id,ivm)
          w%builder_mode = 0
          w%builder_iview = 0
          w%builder_isys = 0
       elseif (win(iview)%vmdata%idx(1) > 0) then
          ! an atom click was delivered: apply the edit and stay in the mode.
          ! Discard stale clicks instead: those delivered while the builder
          ! was not polling (hidden dock tab, collapsed window) or with a
          ! geometry change after the last click-free poll, because the
          ! stored cell index may then refer to a different atom.
          icel = win(iview)%vmdata%idx(1)
          win(iview)%vmdata%idx = 0
          if (glfwGetTime() - w%builder_time < stale_gap .and.&
             sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) &
             call sysc(w%builder_isys)%change_valence(icel,w%builder_mode)
       else
          ! no click pending: stamp the reference time for the stale-click guard
          w%builder_time = glfwGetTime()
       end if
    end if

    ! the valence section
    call iw_text("Valence",highlight=.true.)
    if (iw_button("Increase",disabled=.not.havesys)) &
       call builder_toggle(1)
    call iw_tooltip("Add a hydrogen to each clicked atom in the view,"//&
       " repositioning its terminal substituents (click again to stop)",ttshown)
    if (iw_button("Decrease",sameline=.true.,disabled=.not.havesys)) &
       call builder_toggle(2)
    call iw_tooltip("Remove a hydrogen from each clicked atom in the view,"//&
       " repositioning its terminal substituents (click again to stop)",ttshown)

  contains
    ! Toggle valence mode imode (1 = increase, 2 = decrease): start it on
    ! the main view, or stop it if it is already the active mode.
    subroutine builder_toggle(imode)
      integer, intent(in) :: imode

      integer :: jvm

      jvm = vm_builder_inc
      if (imode == 2) jvm = vm_builder_dec
      if (w%builder_mode == imode) then
         ! stop the active mode and release the view
         call win(w%builder_iview)%viewmode_release_forced(w%id,jvm)
         w%builder_mode = 0
         w%builder_iview = 0
         w%builder_isys = 0
      else
         w%builder_mode = imode
         w%builder_iview = iwin_view
         w%builder_isys = isys
         w%builder_time = glfwGetTime()
         if (imode == 1) then
            call win(iwin_view)%viewmode_set_forced(vm_builder_inc,&
               "Click on atoms to add hydrogens (press any key to exit)...",w%id)
         else
            call win(iwin_view)%viewmode_set_forced(vm_builder_dec,&
               "Click on atoms to remove hydrogens (press any key to exit)...",w%id)
         end if
      end if
    end subroutine builder_toggle
  end subroutine draw_builder

end submodule builder
