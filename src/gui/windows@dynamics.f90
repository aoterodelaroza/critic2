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
    use systems, only: sysc, sys, nsys, sys_init, ok_system, lastchange_geometry
    use dynamics, only: md_dynamics, md_relax
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple, iw_dragfloat_real8,&
       iw_radiobutton, iw_close_event, iw_setpos_bottomright, iw_table_column
    use keybindings, only: get_bind_keyname, BIND_CANCEL
    use tools_io, only: string
    use param, only: kcal2ha, autofs
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, goodparent, ldum, haspress, syschanged
    integer :: isys, iview, isysold
    integer(c_int) :: imode, tflags
    real*8 :: pgpa
    type(ImVec2) :: sz0
    character(len=:), allocatable :: errmsg
    character(len=:,kind=c_char), allocatable, target :: str1

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window acts on the system shown in its anchor view. Keep the outgoing
    ! system: if the view switched, the run on it has to be stopped below
    isysold = w%isys
    goodparent = w%anchor(iview,isys,syschanged)

    ! initialize state
    if (w%firstpass) w%errmsg = ""

    ! resolve the system
    doquit = .not.goodparent
    if (.not.doquit) doquit = .not.associated(win(iview)%sc)
    goodsys = .false.
    if (.not.doquit) goodsys = ok_system(isys,sys_init)

    ! if the view switched systems, stop the run on the one we are leaving
    if (.not.doquit .and. syschanged .and. isysold >= 1 .and. isysold <= nsys) then
       if (ok_system(isysold,sys_init)) then
          if (sysc(isysold)%md_run) call sysc(isysold)%md_stop(errmsg=w%errmsg)
       end if
    end if

    if (goodsys) then
       ! system name
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! method used for MD/relaxation
       call iw_text("Method",highlight=.true.,alignframe=.true.)
       call draw_ff_backend_combo(isys,"##dynamicsengine",21)
       call iw_tooltip("Method for the calculation of energies, forces, and stress.",ttshown)

       ! mode (dynamics vs relaxation), bound live to the run: two radio buttons
       call igAlignTextToFramePadding()
       imode = sysc(isys)%md%mode
       ldum = iw_radiobutton("Dynamics (NVT)",int=imode,intval=int(md_dynamics,c_int))
       call iw_tooltip("Animate the system using an NVT molecular dynamics run",ttshown)
       ldum = iw_radiobutton("Relaxation",int=imode,intval=int(md_relax,c_int),sameline=.true.)
       call iw_tooltip("Relax the system to its nearest energy minimum",ttshown)
       call sysc(isys)%md_set_mode(int(imode))

       ! temperature and timestep
       if (imode == md_dynamics) then
          ! temperature
          ldum = iw_dragfloat_real8("Temperature (K)##dynamicstemp",x1=sysc(isys)%md%temperature,&
             speed=5d0,min=0d0,max=10000d0,decimal=0,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Temperature of the thermostat",ttshown)

          ! speed = timestep (used live)
          ldum = iw_dragfloat_real8("Time step (au)##dynamicsspeed",x1=sysc(isys)%md%dt,&
             speed=0.1d0,min=2d0,max=40d0,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Integration time step for the MD run, in atomic units",ttshown)
       end if

       ! run / pause
       if (.not.sysc(isys)%md_run) then
          if (iw_button("RUN",danger=.true.)) then
             call sysc(isys)%md_start(int(imode),errmsg)
             w%errmsg = errmsg
             if (len_trim(w%errmsg) == 0) &
                win(iview)%forcerender = .true.
          end if
          call iw_tooltip("Start the simulation",ttshown)
       else
          if (iw_button("Stop",danger=.true.)) then
             call sysc(isys)%md_stop(errmsg=w%errmsg)
             win(iview)%forcerender = .true.
          end if
          call iw_tooltip("Stop the simulation ("//&
             trim(get_bind_keyname(BIND_CANCEL))//")",ttshown)
       end if

       ! reset the geometry
       if (iw_button("Reset",sameline=.true.,disabled=.not.sysc(isys)%md%ready)) then
          call sysc(isys)%md%reset(sys(isys)%c)
          sysc(isys)%sc%nextbuildlists_fixcam = .true.
          call sysc(isys)%post_event(lastchange_geometry)
          win(iview)%forcerender = .true.
       end if
       call iw_tooltip("Restore the initial geometry and stop all motion",ttshown)

       ! status: table of live run quantities
       if (sysc(isys)%md%ready) then
          tflags = ImGuiTableFlags_None
          tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
          tflags = ior(tflags,ImGuiTableFlags_RowBg)
          tflags = ior(tflags,ImGuiTableFlags_Borders)
          tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
          str1 = "##dynamicsstatus" // c_null_char
          sz0%x = 0._c_float
          sz0%y = 0._c_float
          if (igBeginTable(c_loc(str1),2,tflags,sz0,0._c_float)) then
             call iw_table_column("Property",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Value",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call igTableHeadersRow()

             ! temperature: MD only (a relaxation has no meaningful temperature)
             if (sysc(isys)%md%mode == md_dynamics) &
                call status_row("Temperature (K)",string(sysc(isys)%md%temperature_now(),'f',decimal=1))
             ! energies: both modes
             call status_row("Energy (Hartree)",string(sysc(isys)%md%epot,'f',decimal=6))
             call status_row("Energy (kcal/mol)",string(sysc(isys)%md%epot/kcal2ha,'f',decimal=3))
             if (sysc(isys)%md%mode == md_relax) then
                ! relaxation: convergence indicator (largest atomic force)
                if (sysc(isys)%md%nat > 0) &
                   call status_row("Max |force| (eV/A)",&
                      string(sysc(isys)%md%maxforce(),'f',decimal=4))
             else
                ! MD: elapsed simulation time and, for crystals, the pressure
                call status_row("Time (fs)",string(sysc(isys)%md%simtime*autofs,'f',decimal=1))
                pgpa = sysc(isys)%md%pressure(sys(isys)%c,haspress)
                if (haspress) &
                   call status_row("Pressure (GPa)",string(pgpa,'f',decimal=3))
             end if

             call igEndTable()
          end if
       end if

       ! errors and auto-stop notices from the run itself
       if (allocated(sysc(isys)%md%errmsg)) then
          if (len_trim(sysc(isys)%md%errmsg) > 0) &
             call iw_text(trim(sysc(isys)%md%errmsg),danger=.true.)
       end if

       ! error message
       if (len_trim(w%errmsg) > 0) &
          call iw_text(trim(w%errmsg),danger=.true.)
    end if

    ! right-align and bottom-align the close button
    call iw_setpos_bottomright(5,1)

    ! close button
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window (the MD run is stopped and freed in window_end)
    if (doquit) call w%end()

  contains
    !> Emit one property/value row of the status table.
    subroutine status_row(prop,val)
      character(len=*), intent(in) :: prop, val
      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
      if (igTableSetColumnIndex(0_c_int)) call iw_text(prop)
      if (igTableSetColumnIndex(1_c_int)) call iw_text(val)
    end subroutine status_row

  end subroutine draw_dynamics

end submodule dynamics
