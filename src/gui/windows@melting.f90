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

! Metal melting demonstration window: melt and refreeze a metal slab or
! nanoparticle with an EAM potential, driving the thermostat temperature live.
submodule (windows) melting
  use interfaces_cimgui
  implicit none

  ! The metals on offer: symbol, atomic number, lattice, conventional
  ! cell parameter (angstrom) and experimental melting point (K).
  integer, parameter :: mt_nmetal = 6
  character(len=2), parameter :: mt_sym(mt_nmetal) = (/"Cu","Au","Al","Ni","Fe","W "/)
  integer, parameter :: mt_z(mt_nmetal) = (/29,79,13,28,26,74/)
  logical, parameter :: mt_bcc(mt_nmetal) = (/.false.,.false.,.false.,.false.,.true.,.true./)
  real*8, parameter :: mt_a(mt_nmetal) = (/3.615d0,4.078d0,4.050d0,3.524d0,2.867d0,3.165d0/)
  real*8, parameter :: mt_tm(mt_nmetal) = (/1358d0,1337d0,933d0,1728d0,1811d0,3695d0/)

  real*8, parameter :: mt_tmax_factor = 1.5d0 ! slider ceiling, in units of the melting point
  real*8, parameter :: mt_vacuum = 15d0 ! vacuum above the slab (angstrom)
  real*8, parameter :: mt_dt = 20d0 ! MD time step (a.u.)
  real*8, parameter :: mt_t0 = 300d0 ! starting temperature (K)
  integer, parameter :: mt_nstep0 = 20 ! starting number of MD steps per frame
  real*8, parameter :: mt_radius_factor = 0.42d0 ! displayed atom radius, in units of the nearest-neighbor distance
  real*8, parameter :: mt_border_factor = 1.2d0 ! atom outline, in units of the usual atom border
  integer, parameter :: mt_nlut = 256 ! colormap look-up table size

  ! constants built on first use: the metal combo options, the order-to-color
  ! look-up table (blue = solid, red = liquid), and the nanoparticle atom count
  ! of the last form settings
  character(len=:,kind=c_char), allocatable :: mt_combo
  logical :: mt_lut_ok = .false.
  real(c_float) :: mt_lut(3,mt_nlut)
  integer :: mt_np_im = 0, mt_np_nshell = 0, mt_np_n = 0

contains

  !> Draw the metal melting demonstration window.
  module subroutine draw_melting(w)
    use systems, only: sysc, sys, sys_init, ok_system, remove_system, lastchange_geometry
    use utils, only: iw_table_headers_row, iw_text, iw_button, iw_tooltip, iw_radiobutton, iw_intstepper,&
       iw_close_event, iw_setpos_bottomright, iw_table_column, iw_combo_simple, iw_dragfloat_real8,&
       file_name_base
    use tools_io, only: string
    use param, only: hartoev, autofs
    class(window), intent(inout), target :: w

    logical :: doquit, goodparent, ldum
    integer :: isys, tflags, im
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
       if (.not.allocated(mt_combo)) then
          mt_combo = ""
          do im = 1, mt_nmetal
             mt_combo = mt_combo // trim(mt_sym(im)) // " ("
             if (mt_bcc(im)) then
                mt_combo = mt_combo // "bcc"
             else
                mt_combo = mt_combo // "fcc"
             end if
             mt_combo = mt_combo // ", melts at " // string(nint(mt_tm(im))) // " K)" // c_null_char
          end do
       end if
       call iw_combo_simple("##mtmetal",mt_combo,w%mt%imetal,sameline=.true.)
       call iw_tooltip("Metal, its crystal structure, and its experimental melting point",ttshown)

       ! geometry
       call iw_text("Geometry",highlight=.true.,alignframe=.true.)
       ldum = iw_radiobutton("Slab",int=w%mt%igeom,intval=0_c_int,sameline=.true.)
       call iw_tooltip("A periodic slab with a free surface on top; the two bottom layers are held &
          &fixed as a crystalline substrate",ttshown)
       ldum = iw_radiobutton("Nanoparticle",int=w%mt%igeom,intval=1_c_int,sameline=.true.)
       call iw_tooltip("A free nanoparticle (cuboctahedron for fcc metals, sphere for bcc metals)",ttshown)

       ! size
       if (w%mt%igeom == 0) then
          ldum = iw_intstepper("mtnx",w%mt%nx,label="Cells (x)",minval=3_c_int,maxval=10_c_int,&
             ndigit=2,tooltip="Number of conventional cells along x")
          ldum = iw_intstepper("mtny",w%mt%ny,label="Cells (y)",minval=3_c_int,maxval=10_c_int,&
             ndigit=2,sameline=.true.,tooltip="Number of conventional cells along y")
          ldum = iw_intstepper("mtnl",w%mt%nlayer,label="Layers",minval=4_c_int,maxval=16_c_int,&
             ndigit=2,sameline=.true.,tooltip="Number of atomic layers in the slab")
       else
          ldum = iw_intstepper("mtns",w%mt%nshell,label="Size",minval=1_c_int,maxval=6_c_int,&
             ndigit=1,tooltip="Size of the nanoparticle: number of atomic shells around the central atom")
          call iw_text("(" // string(mt_np_count(w%mt%imetal+1,int(w%mt%nshell))) // " atoms)",sameline=.true.)
       end if

       ! generate a new system
       if (iw_button("New system",danger=.true.)) then
          if (ok_system(w%isys,sys_init)) call remove_system(w%isys)
          w%isys = 0
          w%errmsg = ""
          call build_melting_system(w%mt%imetal+1,int(w%mt%igeom),int(w%mt%nx),int(w%mt%ny),&
             int(w%mt%nlayer),int(w%mt%nshell),w%isys,w%errmsg)
          ! show the new system in the anchor view (see build_water_cluster)
          if (w%isys > 0) call w%retarget(w%isys)
          ! the system as built: the form may be changed again while it
          ! initializes, and everything below refers to what is on screen
          w%mt%imet = w%mt%imetal + 1
          w%mt%isslab = (w%mt%igeom == 0)
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

             ! temperature
             call iw_text("Temperature",highlight=.true.)
             ldum = iw_dragfloat_real8("##mttemp",x1=sysc(isys)%md%temperature,speed=5d0,min=0d0,&
                max=mt_tmax_factor*mt_tm(im),decimal=0,flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Temperature of the thermostat (K); drag to heat or cool the metal",ttshown)
             call iw_text("K (melts at " // string(nint(mt_tm(im))) // " K)",sameline=.true.)

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
                call status_row("Temperature (K)",string(tnow,'f',decimal=0))
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
      if (w%mt%isslab) then
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
      use representations, only: reptype_bonds, reptype_axes, reptype_atoms, atomborder_def
      integer :: is, i

      is = w%isys
      w%mt%started = .true.
      w%mt%needalign = w%mt%isslab

      ! the run
      sysc(is)%md%temperature = mt_t0
      sysc(is)%md_nstep_frame = mt_nstep0
      call mt_run()

      ! display: no axes and no bonds (the bond list is frozen at the start
      ! of the run and would stretch across the liquid); atoms as large
      ! spheres so the liquid looks dense, outlined so that they can be
      ! told apart where they touch
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
         do i = 1, nat
            w%mt%rgba(1:3,i) = mt_order_color(order(i))
            w%mt%rgba(4,i) = 1._c_float
            if (allocated(sysc(is)%md%frozen)) then
               if (sysc(is)%md%frozen(i)) w%mt%rgba(1:3,i) = rgb_fixed
            end if
         end do
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
      call sysc(is)%sc%show_transient_text(w%id,1,string(nint(tnow)) // " K",rgb_dark,&
         (/0.5d0,0.92d0/),1.5d0)

      ! molten fraction at the bottom, in the color of the atoms at that order
      call sysc(is)%sc%show_transient_text(w%id,2,string(nint(100d0*w%mt%molten)) // "% molten",&
         dark * mt_order_color(1d0-w%mt%molten),(/0.5d0,0.06d0/),1d0)

    end subroutine mt_update_scoreboard

  end subroutine draw_melting

  !xx! private procedures

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
    real*8, allocatable, intent(inout) :: x(:,:)

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

  !> Build the system of the demonstration for metal im: a (100) slab
  !> (igeom = 0) with nx x ny conventional cells in plane and nlayer atomic
  !> layers under a vacuum, periodic; or a nanoparticle (igeom = 1) with
  !> nshell shells, non-periodic. Register it as a new system and return
  !> its id. Errors are returned in errmsg.
  subroutine build_melting_system(im,igeom,nx,ny,nlayer,nshell,id,errmsg)
    use crystalseedmod, only: crystalseed
    use systems, only: add_systems_from_seeds, launch_initialization_thread
    use global, only: rborder_def
    use param, only: isformat_r_derived, bohrtoa
    use tools_io, only: string
    integer, intent(in) :: im, igeom, nx, ny, nlayer, nshell
    integer, intent(inout) :: id
    character(len=:), allocatable, intent(inout) :: errmsg

    type(crystalseed), allocatable :: seed(:)
    integer, allocatable :: idlist(:)
    real*8 :: a, c, xy(2,2)
    integer :: nsite, nat, ix, iy, k, m, i

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

    if (igeom == 0) then
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
    else
       ! nanoparticle, non-periodic (positions in bohr)
       call mt_nanoparticle_sites(mt_bcc(im),nshell,a,seed(1)%x)
       nat = size(seed(1)%x,2)
       seed(1)%useabr = 0
       seed(1)%ismolecule = .true.
       seed(1)%border = rborder_def
       seed(1)%cubic = .false.
       seed(1)%name = trim(mt_sym(im)) // " nanoparticle (" // string(nat) // ")"
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
