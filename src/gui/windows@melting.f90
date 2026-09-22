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

! Metal demonstration window: melt and refreeze a metal slab or nanoparticle,
! or sinter two nanoparticles together, with an EAM potential and the
! thermostat temperature driven live.
submodule (windows) melting
  use interfaces_cimgui
  implicit none

  ! The metals on offer: symbol, atomic number, lattice, conventional
  ! cell parameter (angstrom) and experimental melting point (K).
  integer, parameter :: mt_nmetal = 6
  character(len=2), parameter :: mt_sym(mt_nmetal) = (/"Fe","Cu","Au","Al","Ni","W "/)
  integer, parameter :: mt_z(mt_nmetal) = (/26,29,79,13,28,74/)
  logical, parameter :: mt_bcc(mt_nmetal) = (/.true.,.false.,.false.,.false.,.false.,.true./)
  real*8, parameter :: mt_a(mt_nmetal) = (/2.867d0,3.615d0,4.078d0,4.050d0,3.524d0,3.165d0/)
  real*8, parameter :: mt_tm(mt_nmetal) = (/1811d0,1358d0,1337d0,933d0,1728d0,3695d0/)

  real*8, parameter :: mt_tmax_factor = 2.0d0 ! slider ceiling, in units of the melting point
  real*8, parameter :: mt_vacuum = 15d0 ! vacuum above the slab (angstrom)
  real*8, parameter :: mt_dt = 20d0 ! MD time step (a.u.)
  real*8, parameter :: mt_t0 = 300d0 ! starting temperature (K)
  integer, parameter :: mt_nstep0 = 20 ! starting number of MD steps per frame
  real*8, parameter :: mt_radius_factor = 0.42d0 ! displayed atom radius, in units of the nearest-neighbor distance
  real*8, parameter :: mt_border_factor = 1.2d0 ! atom outline, in units of the usual atom border
  integer, parameter :: mt_nlut = 256 ! colormap look-up table size

  ! The geometries on offer.
  integer, parameter :: mtgeom_slab = 0 ! a periodic slab with a free surface
  integer, parameter :: mtgeom_particle = 1 ! one free nanoparticle
  integer, parameter :: mtgeom_pair = 2 ! two nanoparticles in contact (sintering)

  ! misorientation of the second particle, as the first two Euler angles (degrees)
  real*8, parameter :: mt_pair_rot1 = 37d0
  real*8, parameter :: mt_pair_rot2 = 23d0
  real*8, parameter :: mt_pair_gap = 0.9d0 ! initial surface separation, in nearest-neighbor distances

  ! The temperature units on offer: how each is written, what it is called,
  ! and the affine map from kelvin (x = t * scale + offset).
  integer, parameter :: mt_nunit = 3
  character(len=3), parameter :: mt_uname(0:mt_nunit-1) = &
     (/ character(len=3) :: "K", "°C", "°F" /) ! the degree sign is two bytes in UTF-8
  character(len=18), parameter :: mt_ulong(0:mt_nunit-1) = &
     (/ character(len=18) :: "kelvin", "degrees Celsius", "degrees Fahrenheit" /)
  real*8, parameter :: mt_uscale(0:mt_nunit-1) = (/1d0,1d0,1.8d0/)
  real*8, parameter :: mt_uoffset(0:mt_nunit-1) = (/0d0,-273.15d0,-459.67d0/)
  real(c_float), parameter :: mt_slider_pad = 0.9_c_float ! temperature slider: vertical padding, in text heights
  real(c_float), parameter :: mt_slider_grab = 2.5_c_float ! temperature slider: grab width, in text heights

  ! built on first use: the order-to-color look-up table (blue = solid,
  ! red = liquid) and the nanoparticle atom count of the last form settings
  logical :: mt_lut_ok = .false.
  real(c_float) :: mt_lut(3,mt_nlut)
  integer :: mt_np_im = 0, mt_np_nshell = 0, mt_np_n = 0

contains

  !> Draw the metal melting demonstration window.
  module subroutine draw_melting(w)
    use systems, only: sysc, sys, sys_init, ok_system, remove_system, lastchange_geometry
    use utils, only: iw_table_headers_row, iw_text, iw_button, iw_tooltip, iw_radiobutton, iw_intstepper,&
       iw_close_event, iw_setpos_bottomright, iw_table_column, iw_combo_simple,&
       file_name_base
    use tools_io, only: string
    use param, only: hartoev, autofs
    class(window), intent(inout), target :: w

    logical :: doquit, goodparent, ldum
    integer :: isys, tflags, im, i, ns
    integer(c_int) :: nstep
    real*8 :: eatom, tnow
    integer :: iview
    type(ImVec2) :: sz0
    character(len=:,kind=c_char), allocatable, target :: str1
    character(len=:), allocatable :: errmsg

    logical, save :: ttshown = .false. ! tooltip flag

    ! This window owns the system it generates, so unlike the other view tools
    ! it drives its anchor view rather than following it; all it needs from the
    ! anchor is that it still be a live view.
    iview = w%anchor_view()
    goodparent = (iview > 0)
    doquit = .not.goodparent

    ! initialize state
    if (w%firstpass) w%errmsg = ""

    if (.not.doquit) then
       ! metal
       call iw_text("Metal",highlight=.true.,alignframe=.true.)
       str1 = ""
       do im = 1, mt_nmetal
          str1 = str1 // trim(mt_sym(im)) // " ("
          if (mt_bcc(im)) then
             str1 = str1 // "bcc"
          else
             str1 = str1 // "fcc"
          end if
          str1 = str1 // ", melts at " // mt_tstring(mt_tm(im),w%mt%itempunit) // ")" // c_null_char
       end do
       call iw_combo_simple("##mtmetal",str1,w%mt%imetal,sameline=.true.)
       call iw_tooltip("Metal, its crystal structure, and its experimental melting point",ttshown)

       ! geometry
       call iw_text("Geometry",highlight=.true.,alignframe=.true.)
       ldum = iw_radiobutton("Slab",int=w%mt%igeom,intval=int(mtgeom_slab,c_int),sameline=.true.)
       call iw_tooltip("A periodic slab with a free surface on top; the two bottom layers are held &
          &fixed as a crystalline substrate",ttshown)
       ldum = iw_radiobutton("Nanoparticle",int=w%mt%igeom,intval=int(mtgeom_particle,c_int),&
          sameline=.true.)
       call iw_tooltip("A free nanoparticle (cuboctahedron for fcc metals, sphere for bcc metals)",ttshown)
       ldum = iw_radiobutton("Two particles",int=w%mt%igeom,intval=int(mtgeom_pair,c_int),&
          sameline=.true.)
       call iw_tooltip("Two nanoparticles placed in contact, one of them turned so that the join is a &
          &grain boundary: heat them and watch the neck between them grow (sintering)",ttshown)

       ! size
       if (w%mt%igeom == mtgeom_slab) then
          ldum = iw_intstepper("mtnx",w%mt%nx,label="Cells (x)",minval=3_c_int,maxval=10_c_int,&
             ndigit=2,tooltip="Number of conventional cells along x")
          ldum = iw_intstepper("mtny",w%mt%ny,label="Cells (y)",minval=3_c_int,maxval=10_c_int,&
             ndigit=2,sameline=.true.,tooltip="Number of conventional cells along y")
          ldum = iw_intstepper("mtnl",w%mt%nlayer,label="Layers",minval=4_c_int,maxval=16_c_int,&
             ndigit=2,sameline=.true.,tooltip="Number of atomic layers in the slab")
       elseif (w%mt%igeom == mtgeom_particle) then
          ldum = iw_intstepper("mtns",w%mt%nshell,label="Size",minval=1_c_int,maxval=6_c_int,&
             ndigit=1,tooltip="Size of the nanoparticle: number of atomic shells around the central atom")
          call iw_text("(" // string(mt_np_count(w%mt%imetal+1,int(w%mt%nshell))) // " atoms)",&
             sameline=.true.)
       else
          ! two particles cost twice the atoms, so they are capped smaller to
          ! keep the simulation running at a usable frame rate
          ldum = iw_intstepper("mtnsp",w%mt%nshellp,label="Size",minval=1_c_int,maxval=4_c_int,&
             ndigit=1,tooltip="Size of each nanoparticle: number of atomic shells around its central atom")
          call iw_text("(" // string(2*mt_np_count(w%mt%imetal+1,int(w%mt%nshellp))) // " atoms)",&
             sameline=.true.)
       end if

       ! generate a new system
       if (iw_button("New system",danger=.true.)) then
          if (ok_system(w%isys,sys_init)) call remove_system(w%isys)
          w%isys = 0
          w%errmsg = ""
          ns = int(w%mt%nshell)
          if (w%mt%igeom == mtgeom_pair) ns = int(w%mt%nshellp)
          call build_melting_system(w%mt%imetal+1,int(w%mt%igeom),int(w%mt%nx),int(w%mt%ny),&
             int(w%mt%nlayer),ns,w%isys,w%errmsg)
          ! show the new system in the anchor view (see build_water_cluster)
          if (w%isys > 0) call w%retarget(w%isys)
          ! the system as built: the form may be changed again while it
          ! initializes, and everything below refers to what is on screen
          w%mt%imet = w%mt%imetal + 1
          w%mt%igeom_built = int(w%mt%igeom)
          if (allocated(w%mt%frozen)) deallocate(w%mt%frozen)
          w%mt%started = .false.
       end if
       call iw_tooltip("Build the metal system and start the simulation at room temperature",ttshown)

       ! zoom buttons (touchscreens may not have a mouse wheel)
       if (ok_system(w%isys,sys_init) .and. w%isys == win(iview)%isys) then
          if (iw_button("Zoom +",sameline=.true.)) call mt_zoom(0.15_c_float)
          call iw_tooltip("Zoom in",ttshown)
          if (iw_button("Zoom -",sameline=.true.)) call mt_zoom(-0.15_c_float)
          call iw_tooltip("Zoom out",ttshown)
       end if

       ! auto-start: once the generated system has finished initializing,
       ! set up EAM + NVT dynamics and start the run
       if (ok_system(w%isys,sys_init) .and. .not.w%mt%started) &
          call mt_start()

       ! the running simulation
       isys = w%isys
       if (ok_system(isys,sys_init) .and. w%mt%started) then
          ! a slab is best seen from the side, surface up; the camera can
          ! only be placed once the scene has been built
          if (w%mt%needalign .and. sysc(isys)%sc%isinit >= 2) then
             call sysc(isys)%sc%align_view_axis(2,iup=3)
             win(iview)%forcerender = .true.
             w%mt%needalign = .false.
          end if

          ! run / pause / reset: also when the run is not ready (it stopped
          ! because the structure changed, or a force evaluation failed), so
          ! that Run can initialize it again
          call igSeparator()
          if (sysc(isys)%md_run) then
             if (iw_button("Pause",danger=.true.)) sysc(isys)%md_run = .false.
             call iw_tooltip("Pause the simulation",ttshown)
          else
             if (iw_button("Run",danger=.true.)) call mt_run()
             call iw_tooltip("Resume the simulation",ttshown)
          end if
          if (iw_button("Reset",sameline=.true.,disabled=.not.sysc(isys)%md%ready)) then
             call sysc(isys)%md%reset(sys(isys)%c)
             sysc(isys)%md%temperature = mt_t0
             sysc(isys)%sc%nextbuildlists_fixcam = .true.
             call sysc(isys)%post_event(lastchange_geometry)
             win(iview)%forcerender = .true.
             w%mt%dirty = .true.
          end if
          call iw_tooltip("Restore the initial crystal and go back to room temperature",ttshown)

          if (sysc(isys)%md%ready) then
             im = w%mt%imet

             ! temperature: a slider (not a drag widget) as tall and wide as
             ! the window allows
             call iw_text("Temperature",highlight=.true.)
             call iw_text("(melts at " // mt_tstring(mt_tm(im),w%mt%itempunit) // ")",sameline=.true.)
             call mt_temperature_slider()

             ! units the temperatures are shown in
             call iw_text("Units",highlight=.true.,alignframe=.true.)
             do i = 0, mt_nunit-1
                ldum = iw_radiobutton(trim(mt_uname(i)),int=w%mt%itempunit,intval=int(i,c_int),&
                   sameline=.true.)
                call iw_tooltip("Show the temperatures in " // trim(mt_ulong(i)),ttshown)
             end do

             ! speed
             nstep = int(sysc(isys)%md_nstep_frame,c_int)
             if (iw_intstepper("mtnstep",nstep,label="Steps per frame",minval=1_c_int,&
                maxval=100_c_int,ndigit=3,tooltip="Molecular dynamics steps taken between two &
                &frames: more steps make the simulation run faster, at the cost of a lower frame rate")) &
                sysc(isys)%md_nstep_frame = nstep

             ! local order of every atom, colors, molten fraction
             call mt_color_atoms(isys)
             tnow = sysc(isys)%md%temperature_now()
             eatom = sysc(isys)%md%epot * hartoev / real(sysc(isys)%md%nat,8)
             call mt_update_scoreboard(isys,tnow)

             ! status table
             tflags = ImGuiTableFlags_None
             tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
             tflags = ior(tflags,ImGuiTableFlags_RowBg)
             tflags = ior(tflags,ImGuiTableFlags_Borders)
             tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
             str1 = "##mtstatus" // c_null_char
             sz0%x = 0._c_float
             sz0%y = 0._c_float
             if (igBeginTable(c_loc(str1),2,tflags,sz0,0._c_float)) then
                call iw_table_column("Property",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
                call iw_table_column("Value",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
                call iw_table_headers_row()

                if (allocated(sysc(isys)%md%cl%eam%file)) &
                   call status_row("Potential",file_name_base(sysc(isys)%md%cl%eam%file))
                if (allocated(sysc(isys)%md%frozen)) then
                   call status_row("Atoms (free + fixed)",string(count(.not.sysc(isys)%md%frozen)) //&
                      " + " // string(count(sysc(isys)%md%frozen)))
                else
                   call status_row("Atoms",string(sysc(isys)%md%nat))
                end if
                call status_row("Temperature (" // trim(mt_uname(w%mt%itempunit)) // ")",&
                   string(mt_in_units(tnow,w%mt%itempunit),'f',decimal=0))
                call status_row("Potential energy per atom (eV)",string(eatom,'f',decimal=3))
                call status_row("Molten fraction (%)",string(100d0*w%mt%molten,'f',decimal=0))
                call status_row("Time (ps)",string(sysc(isys)%md%simtime*autofs/1000d0,'f',decimal=2))

                call igEndTable()
             end if

          end if

          ! errors and stop notices from the run itself
          if (allocated(sysc(isys)%md%errmsg)) then
             if (len_trim(sysc(isys)%md%errmsg) > 0) &
                call iw_text(trim(sysc(isys)%md%errmsg),danger=.true.,wrap=.true.)
          end if
       end if

       ! error message
       if (len_trim(w%errmsg) > 0) &
          call iw_text(trim(w%errmsg),danger=.true.,wrap=.true.)
    end if

    ! right-align and bottom-align the close button
    call iw_setpos_bottomright(5,1)

    ! close button
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window (the system and its run are removed in window_end)
    if (doquit) call w%end()

  contains
    !> Emit one property/value row of the status table.
    subroutine status_row(prop,val)
      character(len=*), intent(in) :: prop, val
      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
      if (igTableSetColumnIndex(0_c_int)) call iw_text(prop)
      if (igTableSetColumnIndex(1_c_int)) call iw_text(val)
    end subroutine status_row

    !> The temperature control: a slider spanning the width of the window,
    !> padded to a comfortable height for a finger.
    subroutine mt_temperature_slider()
      use gui_main, only: g, fontsize

      real(c_float) :: t, tmin, tmax
      type(ImVec2) :: pad
      character(len=:,kind=c_char), allocatable, target :: str, sformat

      t = real(mt_in_units(sysc(isys)%md%temperature,w%mt%itempunit),c_float)
      tmin = real(mt_in_units(0d0,w%mt%itempunit),c_float)
      tmax = real(mt_in_units(mt_tmax_factor*mt_tm(im),w%mt%itempunit),c_float)
      str = "##mttemp" // c_null_char
      sformat = "%.0f " // trim(mt_uname(w%mt%itempunit)) // c_null_char

      pad%x = g%Style%FramePadding%x
      pad%y = mt_slider_pad * fontsize%y
      call igPushStyleVar_Vec2(ImGuiStyleVar_FramePadding,pad)
      call igPushStyleVar_Float(ImGuiStyleVar_GrabMinSize,mt_slider_grab * fontsize%y)
      call igPushItemWidth(-1._c_float)
      ! the slider rounds its value to the printed precision, so the bottom
      ! of the bar comes back as a hair below absolute zero (-460 F is
      ! -0.18 K); clamp it, a negative temperature has no meaning
      if (igSliderFloat(c_loc(str),t,tmin,tmax,c_loc(sformat),ImGuiSliderFlags_AlwaysClamp)) &
         sysc(isys)%md%temperature = max(mt_in_kelvin(real(t,8),w%mt%itempunit),0d0)
      call igPopItemWidth()
      call igPopStyleVar(2_c_int)
      call iw_tooltip("Temperature of the thermostat: press anywhere on the bar, or drag it, &
         &to heat or cool the metal",ttshown)

    end subroutine mt_temperature_slider

    !> Zoom the displayed scene in (ratio > 0) or out (ratio < 0) by a fixed
    !> step, for touchscreens with no mouse wheel.
    subroutine mt_zoom(ratio)
      real(c_float), intent(in) :: ratio
      call sysc(w%isys)%sc%cam_zoom(ratio)
      win(iview)%forcerender = .true.
    end subroutine mt_zoom

    !> Start (or resume) the NVT run with the EAM backend. md_start
    !> releases every atom, so the substrate mask goes on after it.
    subroutine mt_run()
      use dynamics, only: md_dynamics
      use energy, only: ff_eam
      use param, only: bohrtoa
      integer :: is, i

      is = w%isys
      sysc(is)%md_backend = ff_eam
      sysc(is)%md_eamfile = "" ! the catalogue default for the metal
      call sysc(is)%md_start(md_dynamics,errmsg)
      if (len_trim(errmsg) > 0) then
         w%errmsg = errmsg
         return
      end if
      ! continuous run
      sysc(is)%md%autostop = .false.
      sysc(is)%md%dt = mt_dt
      ! slab: hold the two bottom layers. md_start releases every atom, so
      ! the mask goes on after it; it is decided once, on the initial
      ! geometry (layer spacing a/2, slab starting at z = 0), because by the
      ! time the run is resumed a liquid atom may have wandered below the cut
      if (w%mt%igeom_built == mtgeom_slab) then
         if (.not.allocated(w%mt%frozen)) then
            allocate(w%mt%frozen(sysc(is)%md%nat))
            do i = 1, sysc(is)%md%nat
               w%mt%frozen(i) = (sysc(is)%md%r0(3,i) < 0.75d0 * mt_a(w%mt%imet) / bohrtoa)
            end do
         end if
         if (size(w%mt%frozen,1) == sysc(is)%md%nat) &
            call sysc(is)%md%set_frozen(w%mt%frozen)
      end if
      w%mt%dirty = .true.
      win(iview)%forcerender = .true.

    end subroutine mt_run

    !> Set up the run on the freshly built system: EAM + NVT at room
    !> temperature, the substrate mask, and the display (no bonds or
    !> axes, atoms as large spheres).
    subroutine mt_start()
      use representations, only: reptype_bonds, reptype_axes, reptype_atoms, atomborder_def,&
         repstyle_ballandstick
      integer :: is, i

      is = w%isys
      w%mt%started = .true.
      w%mt%needalign = (w%mt%igeom_built == mtgeom_slab)

      ! the run
      sysc(is)%md%temperature = mt_t0
      sysc(is)%md_nstep_frame = mt_nstep0
      call mt_run()

      ! a style that draws atoms: above crsmall atoms a new scene defaults to
      ! sticks, which has no atoms object at all, so the view would come up
      ! empty and the per-atom coloring below would have nothing to color
      call sysc(is)%sc%set_style(repstyle_ballandstick)

      ! display: no atoms at the cell edges
      sysc(is)%sc%disp%border = .false.

      ! no axes and no bonds
      sysc(is)%highlight_border = real(mt_border_factor * atomborder_def,c_float)
      do i = 1, sysc(is)%sc%nrep
         if (sysc(is)%sc%rep(i)%type == reptype_axes .or. sysc(is)%sc%rep(i)%type == reptype_bonds) then
            sysc(is)%sc%rep(i)%shown = .false.
         elseif (sysc(is)%sc%rep(i)%type == reptype_atoms) then
            sysc(is)%sc%rep(i)%atoms%radii_type = 2
            sysc(is)%sc%rep(i)%atoms%radii_value = mt_radius_factor * mt_dnn(w%mt%imet)
            call sysc(is)%sc%rep(i)%atoms%style%reset(sysc(is)%sc%rep(i))
         end if
      end do
      sysc(is)%sc%nextbuildlists_fixcam = .true.

    end subroutine mt_start

    !> Color every atom by its local crystalline order (blue = solid,
    !> red = liquid) as a transient highlight, which has to be re-armed
    !> every frame. The order itself is recomputed only when the atoms
    !> move (the run is active) or the state is marked dirty.
    subroutine mt_color_atoms(is)
      use systems, only: atlisttype_ncel_bohr
      integer, intent(in) :: is

      real(c_float), parameter :: rgb_fixed(3) = (/0.55_c_float,0.55_c_float,0.55_c_float/) ! held atoms
      ! the two particles of the sintering geometry: blue and amber stay
      ! apart for the common color-vision deficiencies
      real(c_float), parameter :: rgb_part1(3) = (/0.20_c_float,0.45_c_float,0.75_c_float/)
      real(c_float), parameter :: rgb_part2(3) = (/0.95_c_float,0.65_c_float,0.15_c_float/)
      logical :: bypart

      integer :: nat, i
      integer, allocatable :: idx(:)
      real*8, allocatable :: order(:)
      character(len=:), allocatable :: errmsg

      nat = sysc(is)%md%nat
      if (nat <= 0 .or. nat /= sys(is)%c%ncel) return
      if (allocated(w%mt%rgba)) then
         if (size(w%mt%rgba,2) /= nat) deallocate(w%mt%rgba)
      end if

      if (sysc(is)%md_run .or. w%mt%dirty .or. .not.allocated(w%mt%rgba)) then
         if (.not.allocated(w%mt%rgba)) allocate(w%mt%rgba(4,nat))
         call sysc(is)%md%local_order(sys(is)%c,mt_rcut(w%mt%imet),order,w%mt%molten,errmsg)
         if (len_trim(errmsg) > 0) then
            w%errmsg = errmsg
            return
         end if
         ! the two particles are told apart by color, which is the whole
         ! point of the sintering geometry; everywhere else the color is the
         ! local crystalline order. The builder emits the first particle and
         ! then the second, so the halves of the atom list are the particles
         bypart = (w%mt%igeom_built == mtgeom_pair .and. mod(nat,2) == 0)
         if (bypart) then
            do i = 1, nat
               if (2*i > nat) then
                  w%mt%rgba(1:3,i) = rgb_part2
               else
                  w%mt%rgba(1:3,i) = rgb_part1
               end if
            end do
         else
            do i = 1, nat
               w%mt%rgba(1:3,i) = mt_order_color(order(i))
            end do
         end if
         w%mt%rgba(4,:) = 1._c_float
         if (allocated(sysc(is)%md%frozen)) then
            do i = 1, nat
               if (sysc(is)%md%frozen(i)) w%mt%rgba(1:3,i) = rgb_fixed
            end do
         end if
         w%mt%dirty = .false.
      end if

      idx = (/(i,i=1,nat)/)
      call sysc(is)%highlight_atoms(.true.,idx,atlisttype_ncel_bohr,w%mt%rgba)

    end subroutine mt_color_atoms

    !> Re-arm the on-screen temperature and molten fraction as transient text.
    subroutine mt_update_scoreboard(is,tnow)
      integer, intent(in) :: is
      real*8, intent(in) :: tnow

      real(c_float), parameter :: rgb_dark(3) = (/0.20_c_float,0.20_c_float,0.20_c_float/)
      real(c_float), parameter :: dark = 0.72_c_float ! darken for legibility

      ! temperature at the top
      call sysc(is)%sc%show_transient_text(w%id,1,mt_tstring(tnow,w%mt%itempunit),rgb_dark,&
         (/0.5d0,0.92d0/),1.5d0)

      ! molten fraction at the bottom, in the color of the atoms at that order
      call sysc(is)%sc%show_transient_text(w%id,2,string(nint(100d0*w%mt%molten)) // "% molten",&
         dark * mt_order_color(1d0-w%mt%molten),(/0.5d0,0.06d0/),1d0)

    end subroutine mt_update_scoreboard

  end subroutine draw_melting

  !xx! private procedures

  !> Temperature t (K) converted to the units iunit (0 = kelvin,
  !> 1 = Celsius, 2 = Fahrenheit), for display.
  function mt_in_units(t,iunit) result(x)
    real*8, intent(in) :: t
    integer(c_int), intent(in) :: iunit
    real*8 :: x

    x = t * mt_uscale(iunit) + mt_uoffset(iunit)

  end function mt_in_units

  !> Temperature x, given in the units iunit, back in kelvin.
  function mt_in_kelvin(x,iunit) result(t)
    real*8, intent(in) :: x
    integer(c_int), intent(in) :: iunit
    real*8 :: t

    t = (x - mt_uoffset(iunit)) / mt_uscale(iunit)

  end function mt_in_kelvin

  !> Temperature t (K) written out in the units iunit, rounded, with the
  !> units after it ("1538 °C").
  function mt_tstring(t,iunit) result(str)
    use tools_io, only: string
    real*8, intent(in) :: t
    integer(c_int), intent(in) :: iunit
    character(len=:), allocatable :: str

    str = string(nint(mt_in_units(t,iunit))) // " " // trim(mt_uname(iunit))

  end function mt_tstring

  !> Nearest-neighbor distance (bohr) of metal im.
  function mt_dnn(im) result(dnn)
    use param, only: bohrtoa
    integer, intent(in) :: im
    real*8 :: dnn

    if (mt_bcc(im)) then
       dnn = 0.5d0 * sqrt(3d0) * mt_a(im) / bohrtoa
    else
       dnn = mt_a(im) / (sqrt(2d0) * bohrtoa)
    end if

  end function mt_dnn

  !> Neighbor cutoff (bohr) of the local order parameter for metal im:
  !> the first shell (12) in fcc, the first and second shells (8 + 6) in bcc.
  function mt_rcut(im) result(rcut)
    integer, intent(in) :: im
    real*8 :: rcut

    if (mt_bcc(im)) then
       rcut = 1.3d0 * mt_dnn(im)
    else
       rcut = 1.2d0 * mt_dnn(im)
    end if

  end function mt_rcut

  !> Color of an atom with local order x (0 = liquid, 1 = solid).
  function mt_order_color(x) result(rgb)
    use utils, only: iw_colormap_lut, iw_cmap_rdbu
    real*8, intent(in) :: x
    real(c_float) :: rgb(3)

    integer :: k

    if (.not.mt_lut_ok) then
       call iw_colormap_lut(iw_cmap_rdbu,mt_lut)
       mt_lut_ok = .true.
    end if
    k = 1 + nint(real(mt_nlut-1,8) * min(max(x,0d0),1d0))
    rgb = mt_lut(:,k)

  end function mt_order_color

  !> Number of atoms of the nanoparticle of metal im with nshell shells
  !> (cached: the form is redrawn every frame).
  function mt_np_count(im,nshell) result(n)
    integer, intent(in) :: im, nshell
    integer :: n

    real*8, allocatable :: x(:,:)

    if (im /= mt_np_im .or. nshell /= mt_np_nshell) then
       call mt_nanoparticle_sites(mt_bcc(im),nshell,1d0,x)
       mt_np_im = im
       mt_np_nshell = nshell
       mt_np_n = size(x,2)
    end if
    n = mt_np_n

  end function mt_np_count

  !> Lattice sites of a nanoparticle with nshell shells around a central
  !> atom, in a lattice of cell parameter a (bohr): a cuboctahedron for
  !> fcc (13, 55, 147, ... atoms), a sphere of radius nshell*a for bcc.
  !> Returns the Cartesian coordinates (bohr).
  subroutine mt_nanoparticle_sites(bcc,nshell,a,x)
    logical, intent(in) :: bcc
    integer, intent(in) :: nshell
    real*8, intent(in) :: a
    real*8, allocatable, intent(out) :: x(:,:)

    integer :: i, j, k, n, nmax
    real*8 :: r(3), r2max

    ! sites of both lattices are (i,j,k) * a/2 with: fcc, i+j+k even;
    ! bcc, i, j, k all even or all odd
    if (bcc) then
       nmax = 2*nshell
    else
       nmax = nshell
    end if
    n = 0
    allocate(x(3,(2*nmax+1)**3))
    r2max = (real(nshell,8) * a)**2
    do i = -nmax, nmax
       do j = -nmax, nmax
          do k = -nmax, nmax
             if (bcc) then
                if (mod(abs(i),2) /= mod(abs(j),2) .or. mod(abs(i),2) /= mod(abs(k),2)) cycle
                r = 0.5d0 * a * real((/i,j,k/),8)
                if (dot_product(r,r) > r2max + 1d-8) cycle
             else
                if (mod(abs(i+j+k),2) /= 0) cycle
                ! cuboctahedron: |x|,|y|,|z| <= n a/2 and |x|+|y|+|z| <= n a
                if (abs(i)+abs(j)+abs(k) > 2*nshell) cycle
                r = 0.5d0 * a * real((/i,j,k/),8)
             end if
             n = n + 1
             x(:,n) = r
          end do
       end do
    end do
    x = x(:,1:n)

  end subroutine mt_nanoparticle_sites

  !> Build the system of the demonstration for metal im, one of the
  !> mtgeom_* geometries: a periodic (100) slab with nx x ny conventional
  !> cells in plane and nlayer atomic layers under a vacuum; one
  !> non-periodic nanoparticle of nshell shells; or two such nanoparticles
  !> in contact, the second one turned so that the join is a grain
  !> boundary. Register it as a new system and return its id. Errors are
  !> returned in errmsg.
  subroutine build_melting_system(im,igeom,nx,ny,nlayer,nshell,id,errmsg)
    use crystalseedmod, only: crystalseed
    use systems, only: add_systems_from_seeds, launch_initialization_thread
    use global, only: rborder_def
    use param, only: isformat_r_derived, bohrtoa, rad
    use tools_io, only: string
    use tools_math, only: euler2mat
    integer, intent(in) :: im, igeom, nx, ny, nlayer, nshell
    integer, intent(inout) :: id
    character(len=:), allocatable, intent(inout) :: errmsg

    type(crystalseed), allocatable :: seed(:)
    integer, allocatable :: idlist(:)
    real*8, allocatable :: x1(:,:), y1(:,:)
    real*8 :: a, c, xy(2,2), sep, rot(3,3)
    integer :: nsite, nat, ix, iy, k, m, i, n1

    errmsg = ""
    if (im < 1 .or. im > mt_nmetal) then
       errmsg = "invalid metal"
       return
    end if
    a = mt_a(im) / bohrtoa

    allocate(seed(1))
    seed(1)%nspc = 1
    allocate(seed(1)%spc(1))
    seed(1)%spc(1)%z = mt_z(im)
    seed(1)%spc(1)%name = trim(mt_sym(im))

    if (igeom == mtgeom_slab) then
       ! (100) slab: layers a/2 apart; each layer is a square lattice with
       ! the conventional cell, and consecutive layers alternate between
       ! the two motifs of the conventional cell
       if (mt_bcc(im)) then
          nsite = 1
       else
          nsite = 2
       end if
       nat = nsite * nx * ny * nlayer
       c = 0.5d0 * a * real(nlayer,8) + mt_vacuum / bohrtoa
       allocate(seed(1)%x(3,nat))
       i = 0
       do k = 0, nlayer-1
          ! in-plane motif of this layer, in units of the conventional cell
          if (mt_bcc(im)) then
             if (mod(k,2) == 0) then
                xy(:,1) = (/0d0,0d0/)
             else
                xy(:,1) = (/0.5d0,0.5d0/)
             end if
          else
             if (mod(k,2) == 0) then
                xy(:,1) = (/0d0,0d0/)
                xy(:,2) = (/0.5d0,0.5d0/)
             else
                xy(:,1) = (/0.5d0,0d0/)
                xy(:,2) = (/0d0,0.5d0/)
             end if
          end if
          do ix = 0, nx-1
             do iy = 0, ny-1
                do m = 1, nsite
                   i = i + 1
                   seed(1)%x(1,i) = (real(ix,8) + xy(1,m)) / real(nx,8)
                   seed(1)%x(2,i) = (real(iy,8) + xy(2,m)) / real(ny,8)
                   seed(1)%x(3,i) = 0.5d0 * a * real(k,8) / c
                end do
             end do
          end do
       end do
       seed(1)%m_x2c = 0d0
       seed(1)%m_x2c(1,1) = a * real(nx,8)
       seed(1)%m_x2c(2,2) = a * real(ny,8)
       seed(1)%m_x2c(3,3) = c
       seed(1)%useabr = 2
       seed(1)%ismolecule = .false.
       seed(1)%name = trim(mt_sym(im)) // "(100) slab " // string(nx) // "x" // string(ny) //&
          "x" // string(nlayer)
    elseif (igeom == mtgeom_particle) then
       ! nanoparticle, non-periodic (positions in bohr)
       call mt_nanoparticle_sites(mt_bcc(im),nshell,a,seed(1)%x)
       nat = size(seed(1)%x,2)
       seed(1)%cubic = .false.
       seed(1)%name = trim(mt_sym(im)) // " nanoparticle (" // string(nat) // ")"
    elseif (igeom == mtgeom_pair) then
       ! two nanoparticles in contact along x, non-periodic (positions in
       ! bohr). The second one is turned to a generic orientation: two copies
       ! in register would simply be one crystal, and the join has to be a
       ! grain boundary for there to be any sintering to watch. A turn about
       ! the contact axis alone would not do, since it leaves the facing
       ! facets parallel and the contact flat
       call mt_nanoparticle_sites(mt_bcc(im),nshell,a,x1)
       n1 = size(x1,2)
       nat = 2 * n1
       rot = euler2mat((/mt_pair_rot1 * rad, mt_pair_rot2 * rad, 0d0/))
       y1 = matmul(rot,x1)
       ! the gap is measured on the coordinates themselves: the turned copy
       ! reaches further along x than the original, so any radius formula
       ! would push the two into each other
       sep = mt_pair_gap * mt_dnn(im) + maxval(x1(1,:)) - minval(y1(1,:))
       allocate(seed(1)%x(3,nat))
       do i = 1, n1
          seed(1)%x(:,i) = x1(:,i) - (/0.5d0*sep,0d0,0d0/)
          seed(1)%x(:,n1+i) = y1(:,i) + (/0.5d0*sep,0d0,0d0/)
       end do
       ! a cube: the pair's own box is long and thin, and would be re-fitted
       ! (rebuilding the cell and the environment) as the pair tumbles
       seed(1)%cubic = .true.
       seed(1)%name = trim(mt_sym(im)) // " two nanoparticles (2x" // string(n1) // ")"
    else
       errmsg = "invalid geometry"
       return
    end if

    ! the two free-standing geometries share their seed conventions
    if (igeom /= mtgeom_slab) then
       seed(1)%useabr = 0
       seed(1)%ismolecule = .true.
       seed(1)%border = rborder_def
    end if

    ! common: one species, built in memory (no source file), no symmetry
    seed(1)%nat = nat
    allocate(seed(1)%is(nat),seed(1)%atname(nat))
    seed(1)%is = 1
    seed(1)%atname = trim(mt_sym(im))
    seed(1)%havesym = 0
    seed(1)%findsym = 0
    seed(1)%neqv = 0
    seed(1)%ncv = 0
    seed(1)%havex0 = .false.
    seed(1)%molx0 = 0d0
    seed(1)%isused = .true.
    seed(1)%isformat = isformat_r_derived
    seed(1)%file = ""

    ! register the system; add_systems_from_seeds selects it in the tree, so the
    ! main view displays it
    call add_systems_from_seeds(1,seed,idlist=idlist)
    call launch_initialization_thread()
    id = 0
    if (allocated(idlist)) then
       if (size(idlist) >= 1) id = idlist(1)
    end if

  end subroutine build_melting_system

end submodule melting
