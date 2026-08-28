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

! Routines for the rattle structure window: generate new systems from the
! current one by displacing its atoms randomly or by sampling a molecular
! dynamics run. This is the GUI counterpart of WRITE BULK RATTLE and
! WRITE BULK MD; the random displacements share their generator with the
! keyword (bulk_rattle_seeds in crystalmod).
!
! The dynamics is not run in one blocking call (bulk_md_seeds in dynamics,
! which the keyword uses). Instead the window drives the parent system's own
! interactive run: the main loop advances it one step per frame, so the
! sampling animates in the view and the interface stays responsive.
submodule (windows) rattle
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! how the new structures are generated
  integer(c_int), parameter :: rm_random = 0 ! random displacements (RATTLE)
  integer(c_int), parameter :: rm_md = 1 ! snapshots from an MD run (MD)

  ! settings remembered from the last successful run (this session)
  integer(c_int) :: lastmode = rm_random
  integer(c_int) :: lastnstruct = 10_c_int
  real(c_float) :: lastmag = 0.02_c_float
  integer(c_int) :: lastnini = 1000_c_int
  integer(c_int) :: lastngen = 100_c_int
  integer(c_int) :: lastnstride = 1_c_int
  logical :: lastcloseafter = .true.

contains

  !> Draw the rattle structure window: generate a group of new systems from
  !> the system in the anchor view, either by rattling the atoms at random or
  !> by collecting snapshots from a molecular-dynamics run.
  module subroutine draw_rattle(w)
    use systems, only: sysc, sys, ok_system, sys_init, add_systems_from_seeds,&
       add_group_master, launch_initialization_thread
    use crystalmod, only: bulk_rattle_seeds
    use dynamics, only: md_dynamics, md_benchmark
    use utils, only: iw_text, iw_button, iw_tooltip, iw_radiobutton, iw_checkbox,&
       iw_close_event, iw_setpos_bottomright, iw_intstepper, iw_dragfloat_realc,&
       iw_dragfloat_real8, duration_string
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_OK_FOCUSED_DIALOG, BIND_CANCEL
    use tools_io, only: string, uout
    use param, only: autofs
    class(window), intent(inout), target :: w

    logical :: doquit, ldum, ok, syschanged, running
    integer :: isys, iview, nsteptot
    character(len=:), allocatable :: errmsg

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    if (.not.doquit) doquit = .not.ok_system(isys,sys_init)

    ! initialize state
    if (w%firstpass) then
       w%rattle_mode = lastmode
       w%rattle_nstruct = lastnstruct
       w%rattle_mag = lastmag
       w%rattle_nini = lastnini
       w%rattle_ngen = lastngen
       w%rattle_nstride = lastnstride
       w%rattle_closeafter = lastcloseafter
       w%rattle_running = .false.
    end if
    if (w%firstpass .or. syschanged) then
       call clear_msgs()
       w%rattle_tinit = -1d0
    end if

    ! the cancel keybinding ends a sampling run, keeping what it collected.
    ! When the view is the focused target the global handler has already
    ! stopped the dynamics by now, and advance_run() picks that up instead
    if (w%rattle_running) then
       if (is_bind_event(BIND_CANCEL,norepeat=.true.)) call finish_run(.true.)
    end if

    ! advance a sampling run in progress: collect the snapshot that came due
    ! this frame, and close the run when the last one is in
    if (w%rattle_running) call advance_run()
    running = w%rattle_running

    ! the system the structures are generated from
    if (.not.doquit) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)
    end if

    ! how the structures are made; the outcome of a run does not carry over
    ! to the other method, whose parameters it says nothing about. The
    ! parameters are fixed while a run is going
    call igBeginDisabled(logical(running,c_bool))
    call iw_text("Method",highlight=.true.)
    if (iw_radiobutton("Random displacements##rattlemode",int=w%rattle_mode,intval=rm_random)) &
       call clear_msgs()
    call iw_tooltip("Displace every atom of the structure by the same distance in a random&
       & direction",ttshown)
    if (iw_radiobutton("Molecular dynamics##rattlemode",int=w%rattle_mode,intval=rm_md,&
       sameline=.true.)) call clear_msgs()
    call iw_tooltip("Run an NVT constant-temperature molecular dynamics on this system and extract&
       & snapshots from it",ttshown)

    ! the parameters of whichever method is selected
    call iw_text("Details",highlight=.true.)
    if (w%rattle_mode == rm_random) then
       call iw_text("Structures",alignframe=.true.)
       ldum = iw_intstepper("nstruct##rattlenstruct",w%rattle_nstruct,minval=1_c_int,&
          maxval=10000_c_int,notlive=.true.,sameline=.true.,&
          tooltip="Number of rattled structures to generate")

       call iw_text("Displacement",alignframe=.true.)
       ldum = iw_dragfloat_realc("(Å)##rattlemag",x1=w%rattle_mag,speed=0.005_c_float,&
          min=0.0001_c_float,max=10._c_float,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp,&
          sameline=.true.)
       call iw_tooltip("Distance every atom is moved; the direction is random and&
          & different for each atom",ttshown)
    elseif (.not.doquit) then
       call iw_text("Force field",alignframe=.true.)
       call draw_ff_backend_combo(isys,"##rattleff",21)

       ! the run parameters live on the system, shared with the dynamics window
       call iw_text("Temperature",alignframe=.true.)
       ldum = iw_dragfloat_real8("(K)##rattletemp",x1=sysc(isys)%md%temperature,&
          speed=1d0,min=0d0,max=1d4,decimal=2,&
          flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
       call iw_tooltip("Temperature of the NVT dynamics",ttshown)

       ! shown in fs, stored in atomic units (scale), which is what the
       ! dynamics and the STEP option to WRITE BULK MD use
       call iw_text("Time step",alignframe=.true.)
       ldum = iw_dragfloat_real8("(fs)##rattledt",x1=sysc(isys)%md%dt,speed=0.005d0,&
          min=0.01d0,max=1d0,scale=autofs,decimal=4,&
          flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
       call iw_tooltip("Integration time step of the dynamics",ttshown)
       call iw_text("(= " // string(sysc(isys)%md%dt,'f',decimal=2) // " a.u.)",sameline=.true.)

       call iw_text("Equilibration",alignframe=.true.)
       ldum = iw_intstepper("nini##rattlenini",w%rattle_nini,minval=0_c_int,&
          maxval=1000000_c_int,notlive=.true.,sameline=.true.,&
          tooltip="Number of equilibration steps run before the first snapshot is taken")

       call iw_text("Snapshots",alignframe=.true.)
       ldum = iw_intstepper("ngen##rattlengen",w%rattle_ngen,minval=1_c_int,&
          maxval=10000_c_int,notlive=.true.,sameline=.true.,&
          tooltip="Number of snapshots collected, one every stride steps")

       call iw_text("Stride",alignframe=.true.)
       ldum = iw_intstepper("nstride##rattlenstride",w%rattle_nstride,minval=1_c_int,&
          maxval=100000_c_int,notlive=.true.,sameline=.true.,&
          tooltip="Number of dynamics steps between two consecutive snapshots")

       ! the cost of the run, which is what the user is really choosing: the
       ! dynamics advances one step per frame, so this is also a frame count
       nsteptot = int(w%rattle_nini) + int(w%rattle_ngen) * int(w%rattle_nstride)
       call iw_text("(" // string(nsteptot) // " steps in total)",alignframe=.true.)

       ! measure what a step of this dynamics costs, and turn it into an
       ! estimate for the whole run. The measurement is kept, so the estimate
       ! follows the step count as the user changes it
       if (w%rattle_estbackend /= sysc(isys)%md_backend) w%rattle_tinit = -1d0
       if (iw_button("Estimate cost##rattleestimate",sameline=.true.)) call estimate_cost()
       call iw_tooltip("Run a few steps of this dynamics and discard the results to&
          & measure the cost and estimate how long the whole run will take",ttshown)
       if (w%rattle_tinit >= 0d0) then
          call iw_text("Estimated cost: " //&
             duration_string(w%rattle_tinit + nsteptot * w%rattle_tstep) // " (" //&
             duration_string(w%rattle_tstep) // " per step)",highlight=.true.)
       end if
    end if
    call igEndDisabled()

    ! progress of a run in progress, in its own line
    if (running .and. .not.doquit) then
       if (sysc(isys)%md%istep - w%rattle_istep0 < int(w%rattle_nini)) then
          call iw_text("Equilibrating: step " //&
             string(sysc(isys)%md%istep - w%rattle_istep0) // " of " //&
             string(w%rattle_nini),highlight=.true.)
       else
          call iw_text("Sampling: " // string(w%rattle_ncoll) // " of " //&
             string(w%rattle_ngen) // " snapshots",highlight=.true.)
       end if
    end if

    ! the outcome of the last run
    if (len_trim(w%errmsg) > 0) then
       call iw_text(w%errmsg,danger=.true.,wrap=.true.)
    elseif (len_trim(w%okmsg) > 0) then
       call iw_text(w%okmsg,highlight=.true.,wrap=.true.)
    end if

    ! right-align and bottom-align for the rest of the contents: the checkbox
    ! plus the two buttons of the row below, whose first one is Stop while a
    ! run is going and Generate otherwise
    call iw_setpos_bottomright(merge(31,35,running),2,ncheck=1)

    ! close after generating, next to the button it applies to
    call igBeginDisabled(logical(running,c_bool))
    ldum = iw_checkbox("Close after generating##rattlecloseafter",w%rattle_closeafter)
    call iw_tooltip("If checked, generating the structures closes this window",ttshown)
    call igEndDisabled()

    ! final buttons: generate, or stop a run in progress
    if (running) then
       if (iw_button("Stop",sameline=.true.,danger=.true.)) call finish_run(.true.)
       call iw_tooltip("Stop the dynamics now and keep the snapshots collected so far ("//&
          trim(get_bind_keyname(BIND_CANCEL))//")",ttshown)
    else
       ok = iw_button("Generate",disabled=doquit,sameline=.true.)
       ok = ok .or. (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG) .and. .not.doquit)
       if (ok) then
          call clear_msgs()
          if (w%rattle_mode == rm_random) then
             call generate_rattle()
          else
             call start_run()
          end if
       end if
       call iw_tooltip("Generate the new structures and add them to the tree as a&
          & group of new systems",ttshown)
    end if

    ! final buttons: close
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  contains
    ! Forget the outcome of the previous run.
    subroutine clear_msgs()

      w%errmsg = ""
      w%okmsg = ""

    end subroutine clear_msgs

    ! Measure the cost of one step of the dynamics
    subroutine estimate_cost()
      character(len=:), allocatable :: errmsg
      integer :: nstep
      real*8 :: tinit, tstep

      call clear_msgs()
      w%rattle_tinit = -1d0
      call md_benchmark(sys(isys)%c,sysc(isys)%md_backend,md_dynamics,&
         sysc(isys)%md%temperature,sysc(isys)%md%dt,5,0.5d0,nstep,tinit,tstep,errmsg)
      if (len_trim(errmsg) > 0) then
         w%errmsg = errmsg
         return
      end if
      w%rattle_tinit = tinit
      w%rattle_tstep = tstep
      w%rattle_estbackend = sysc(isys)%md_backend

    end subroutine estimate_cost

    ! Remember the settings of a run that completed, and close the window if
    ! that is what the user asked for.
    subroutine run_done()

      lastmode = w%rattle_mode
      lastnstruct = w%rattle_nstruct
      lastmag = w%rattle_mag
      lastnini = w%rattle_nini
      lastngen = w%rattle_ngen
      lastnstride = w%rattle_nstride
      lastcloseafter = w%rattle_closeafter
      if (w%rattle_closeafter) doquit = .true.

    end subroutine run_done

    ! Add the seeds in seed(1:nseed) to the tree as a group of their own,
    ! headed by an entry labelled label.
    ! is = the system the structures were made from. It is not always the one
    ! the anchor view shows now: a sampling run keeps its own, and the view may
    ! have moved on or gone away (isys = 0) by the time the run ends
    subroutine add_group(is,label,seed,nseed)
      integer, intent(in) :: is
      character(len=*), intent(in) :: label
      type(crystalseed), allocatable, intent(in) :: seed(:)
      integer, intent(in) :: nseed

      integer :: imas
      character(len=:), allocatable :: parent

      parent = "unknown"
      if (ok_system(is,sys_init)) parent = trim(sysc(is)%seed%name)

      call add_group_master(label,parent,imas)
      call add_systems_from_seeds(nseed,seed,imaster=imas,noselect=.true.)
      call launch_initialization_thread()

      write (uout,'("Generated ",A," structures from system ",A,": ",A," (",A,")")') &
         string(nseed), string(is), parent, label
      write (uout,*)

    end subroutine add_group

    ! Generate the randomly displaced structures. This is immediate: no
    ! energy is evaluated, so there is nothing to watch.
    subroutine generate_rattle()
      type(crystalseed), allocatable :: seed(:)
      integer :: nseed, i

      nseed = 0
      allocate(seed(w%rattle_nstruct))
      call bulk_rattle_seeds(sys(isys)%c,int(w%rattle_nstruct),&
         real(w%rattle_mag,8)/bohrtoa,seed,nseed)
      if (nseed == 0) then
         w%errmsg = "No structures were generated"
         return
      end if

      do i = 1, nseed
         seed(i)%name = trim(sysc(isys)%seed%name) // " (rattled " // string(i) // ")"
      end do
      call add_group(isys,"rattled, " // string(real(w%rattle_mag,8),'f',decimal=3) // " Å",&
         seed,nseed)
      w%okmsg = "Generated " // string(nseed) // " rattled structures"
      call run_done()

    end subroutine generate_rattle

    ! Start the sampling dynamics on the parent system. From here on the main
    ! loop advances it one step per frame and advance_run() collects.
    subroutine start_run()

      ! the run is the system's own, and the dynamics window may already be
      ! using it; sampling would stop it and restore the geometry underneath
      if (sysc(isys)%md_run) then
         w%errmsg = "A dynamics run is already active on this system"
         return
      end if
      call sysc(isys)%md_start(md_dynamics,errmsg)
      if (len_trim(errmsg) > 0) then
         w%errmsg = errmsg
         return
      end if
      ! md_start reuses a run that is still ready, whose restore point is the
      ! geometry *it* started from; rebase it on the structure as it is now,
      ! so finishing the sampling puts this geometry back
      call sysc(isys)%md%rebase(sys(isys)%c)

      if (allocated(w%rattle_seed)) deallocate(w%rattle_seed)
      allocate(w%rattle_seed(w%rattle_ngen))
      w%rattle_running = .true.
      w%rattle_isys = isys
      w%rattle_istep0 = sysc(isys)%md%istep
      w%rattle_ncoll = 0
      w%rattle_nextdue = int(w%rattle_nini) + int(w%rattle_nstride)
      w%rattle_runtemp = sysc(isys)%md%temperature

    end subroutine start_run

    ! Collect the snapshot that came due this frame, and end the run when the
    ! last one is in. Also ends it if the dynamics stopped underneath us: the
    ! run is the system's own, so a force-field failure, a structure edit or
    ! the global cancel keybinding can all stop it.
    subroutine advance_run()
      integer :: nstep

      ! the view moved to another system: this run is no longer the one on
      ! screen, so end it (closing the window is handled by window_end)
      if (w%rattle_isys /= isys) then
         call finish_run(.true.)
         return
      end if
      if (.not.sysc(isys)%md_run) then
         ! the run is no longer going: either something stopped it from outside
         ! (the cancel keybinding, another window), which is not an error, or
         ! the dynamics itself failed, which is
         if (allocated(sysc(isys)%md%errmsg)) then
            if (len_trim(sysc(isys)%md%errmsg) > 0) &
               w%errmsg = sysc(isys)%md%fail_message("The dynamics stopped before the run was complete")
         end if
         call finish_run(.true.)
         return
      end if

      ! the dynamics steps every frame but this routine only runs while the
      ! window is drawn, so the due step can be overshot (a collapsed window,
      ! a docked tab in the background). Track the next due step rather than
      ! testing for an exact multiple, and re-anchor it on the step actually
      ! taken, so consecutive snapshots stay at least nstride steps apart
      nstep = sysc(isys)%md%istep - w%rattle_istep0
      if (nstep >= w%rattle_nextdue .and. w%rattle_ncoll < int(w%rattle_ngen)) then
         w%rattle_ncoll = w%rattle_ncoll + 1
         w%rattle_nextdue = nstep + int(w%rattle_nstride)
         call sys(isys)%c%makeseed(w%rattle_seed(w%rattle_ncoll),.false.)
         w%rattle_seed(w%rattle_ncoll)%name = trim(sysc(w%rattle_isys)%seed%name) //&
            " (MD step " // string(nstep) // ")"
      end if
      if (w%rattle_ncoll >= int(w%rattle_ngen)) call finish_run(.false.)

    end subroutine advance_run

    ! End the sampling run: put the structure back where it started, stop the
    ! dynamics, and add whatever was collected to the tree. stopped = the run
    ! was cut short (by the user or by the window closing).
    subroutine finish_run(stopped)
      logical, intent(in) :: stopped

      integer :: is

      is = w%rattle_isys
      w%rattle_running = .false.
      if (ok_system(is,sys_init)) then
         ! the sampling must not move the structure it sampled from
         if (sysc(is)%md%ready) call sysc(is)%md%reset(sys(is)%c)
         call sysc(is)%md_stop(errmsg)
         if (len_trim(errmsg) > 0) w%errmsg = errmsg
      end if

      if (w%rattle_ncoll > 0) then
         call add_group(is,"NVT MD run, " // string(nint(w%rattle_runtemp)) // " K",&
            w%rattle_seed,w%rattle_ncoll)
         if (stopped) then
            w%okmsg = "Stopped after " // string(w%rattle_ncoll) // " snapshots"
         else
            w%okmsg = "Generated " // string(w%rattle_ncoll) // " snapshots"
         end if
      elseif (stopped .and. len_trim(w%errmsg) == 0) then
         w%okmsg = "Stopped before the first snapshot"
      end if
      if (allocated(w%rattle_seed)) deallocate(w%rattle_seed)
      if (.not.stopped .and. len_trim(w%errmsg) == 0) call run_done()

    end subroutine finish_run

  end subroutine draw_rattle

end submodule rattle
