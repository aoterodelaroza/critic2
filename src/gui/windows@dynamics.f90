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

! Routines for the interactive molecular-dynamics window.
submodule (windows) dynamics
  use interfaces_cimgui
  implicit none
contains

  !> Draw the interactive dynamics window: run a real-time MD or relaxation of
  !> the parent view's system, with controls for temperature, speed, engine, and
  !> the ability to grab and drag atoms in the view.
  module subroutine draw_dynamics(w)
    use scenes, only: scene
    use systems, only: sysc, sys, sys_init, ok_system, lastchange_geometry
    use energy, only: ff_builtin, ff_tblite
    use dynamics, only: md_dynamics, md_relax
    use utils, only: iw_text, iw_button, iw_tooltip, iw_calcwidth, iw_combo_simple,&
       iw_dragfloat_real8
    use gui_main, only: g
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use tools_io, only: string
    class(window), intent(inout), target :: w

    type(scene), pointer :: sc
    logical :: doquit, goodsys, goodparent, needinit, ldum
    integer :: isys, imode, ibackend
    type(ImVec2) :: szavail
    character(len=:), allocatable :: errmsg

    character(len=*), parameter :: str_backend = &
       "Built-in force field" // c_null_char // "GFN2-xTB (tblite)" // c_null_char
    character(len=*), parameter :: str_mode = &
       "Temperature (MD)" // c_null_char // "Relaxation" // c_null_char

    logical, save :: ttshown = .false. ! tooltip flag

    ! do we have a good parent window (a view)?
    goodparent = w%idparent > 0 .and. w%idparent <= nwin
    if (goodparent) goodparent = win(w%idparent)%isinit
    if (goodparent) goodparent = (win(w%idparent)%type == wintype_view)

    ! initialize state
    if (w%firstpass) w%errmsg = ""

    ! resolve the system and its scene
    isys = 0
    doquit = .not.goodparent
    nullify(sc)
    if (.not.doquit) then
       if (associated(win(w%idparent)%sc)) then
          sc => win(w%idparent)%sc
          isys = sc%id
       else
          doquit = .true.
       end if
    end if
    goodsys = .false.
    if (.not.doquit) goodsys = ok_system(isys,sys_init)

    if (goodsys) then
       ! system name
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! energy engine
       call igAlignTextToFramePadding()
       call iw_text("Engine",highlight=.true.)
       ibackend = sc%md_backend
       call igSameLine(0._c_float,-1._c_float)
       call igPushItemWidth(iw_calcwidth(21,1))
       call iw_combo_simple("##dynamicsengine",str_backend,ibackend)
       call igPopItemWidth()
       sc%md_backend = ibackend
       call iw_tooltip("Energy and force engine. The built-in force field always works; GFN2-xTB &
          &requires critic2 compiled with tblite support.",ttshown)

       ! mode (dynamics vs relaxation), bound live to the run
       call igAlignTextToFramePadding()
       call iw_text("Mode",highlight=.true.)
       imode = sc%md%mode
       call igSameLine(0._c_float,-1._c_float)
       call igPushItemWidth(iw_calcwidth(17,1))
       call iw_combo_simple("##dynamicsmode",str_mode,imode)
       call igPopItemWidth()
       sc%md%mode = imode
       call iw_tooltip("Temperature: animate the system with a thermostat at the chosen temperature. &
          &Relaxation: drive the system to its nearest energy minimum.",ttshown)

       ! temperature (used live by the thermostat)
       ldum = iw_dragfloat_real8("Temperature (K)##dynamicstemp",x1=sc%md%temperature,&
          speed=1d0,min=0d0,max=2000d0,decimal=0,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Target temperature of the thermostat, in kelvin",ttshown)

       ! speed = timestep (used live)
       ldum = iw_dragfloat_real8("Speed##dynamicsspeed",x1=sc%md%dt,&
          speed=0.1d0,min=2d0,max=40d0,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Simulation speed (integration timestep, atomic units)",ttshown)

       ! run / pause
       if (.not.sc%md_run) then
          if (iw_button("Run")) then
             needinit = (.not.sc%md%ready) .or. (sc%md%cl%backend /= sc%md_backend)
             if (needinit) then
                if (allocated(errmsg)) deallocate(errmsg)
                call sc%md%init(sys(isys)%c,backend=sc%md_backend,errmsg=errmsg)
                sc%md%mode = imode
             end if
             if (allocated(errmsg)) then
                w%errmsg = errmsg
             else
                w%errmsg = ""
                sc%md_run = .true.
                win(w%idparent)%forcerender = .true.
             end if
          end if
          call iw_tooltip("Start (or resume) the simulation",ttshown)
       else
          if (iw_button("Pause",danger=.true.)) then
             sc%md_run = .false.
             win(w%idparent)%forcerender = .true.
          end if
          call iw_tooltip("Pause the simulation",ttshown)
       end if

       ! reset the geometry
       if (iw_button("Reset",sameline=.true.,disabled=.not.sc%md%ready)) then
          call sc%md%reset(sys(isys)%c)
          sc%nextbuildlists_fixcam = .true.
          call sysc(isys)%post_event(lastchange_geometry)
          win(w%idparent)%forcerender = .true.
       end if
       call iw_tooltip("Restore the initial geometry and stop all motion",ttshown)

       ! status line
       if (sc%md%ready) then
          call iw_text("T = " // string(sc%md%temperature_now(),'f',decimal=0) // " K" //&
             "   E = " // string(sc%md%epot,'f',decimal=5) // " Ha")
       end if

       ! usage hint
       call iw_text("Drag an atom with the left mouse button to pull it around.")

       ! error message
       if (len_trim(w%errmsg) > 0) &
          call iw_text(trim(w%errmsg),danger=.true.)
    end if

    ! right-align and bottom-align the close button
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(5,1,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! close button
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit = close the window (the MD run is stopped and freed in window_end)
    if (doquit) call w%end()

  end subroutine draw_dynamics

end submodule dynamics
