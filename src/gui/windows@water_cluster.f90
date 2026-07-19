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

! Water-cluster demonstration window: an interactive TIP4P relaxation game.
! Generates a cluster of N water molecules (random or in a row), FIRE-relaxes
! it with the built-in TIP4P model, and shows the binding energy on screen
! recolored by its closeness to the known global minimum for that cluster size.
submodule (windows) water_cluster
  use interfaces_cimgui
  implicit none

  ! Global-minimum energies (kJ/mol), indexed by the number of water
  ! molecules. These are critic2's own TIP4P energies at each
  ! Cambridge Cluster Database global-minimum geometry (D. J. Wales
  ! and M. P. Hodges, Chem. Phys.  Lett. 286, 65 (1998)). critic2's
  ! TIP4P adds harmonic q-TIP4P/F intramolecular restraints
  ! (energy@proc.F90:45-58), so the monomers flex and bind ~4-5% more
  ! than the rigid TIP4P Wales tabulated.
  integer, parameter :: wc_nmin = 2        ! smallest cluster with a reference
  integer, parameter :: wc_nmax = 21       ! largest cluster with a reference
  integer, parameter :: wc_ngen_max = 100  ! largest cluster the user can generate
  real*8, parameter :: wc_record_kjmol(wc_nmin:wc_nmax) = (/ &
     -26.84746d0,  -72.94662d0, -122.22162d0, -159.44335d0, -207.12239d0, &
     -255.29024d0, -320.18069d0, -360.80033d0, -409.44119d0, -451.80108d0, &
     -516.69115d0, -558.32278d0, -610.71658d0, -658.05751d0, -713.99910d0, &
     -758.57272d0, -810.16291d0, -860.34569d0, -914.82736d0, -961.23746d0/)

  ! TIP4P reference monomer geometry (must match src/energy@proc.F90).
  real*8, parameter :: wc_roh_ang = 0.9572d0 ! reference O-H distance (angstrom)
  real*8, parameter :: wc_hoh_deg = 104.52d0 ! reference H-O-H angle (degrees)

contains

  !> Draw the water-cluster demonstration window: controls to generate a
  !> cluster of N waters (random or in a row), a continuously running TIP4P
  !> FIRE relaxation, and an on-screen binding-energy scoreboard colored by
  !> closeness to the global minimum.
  module subroutine draw_water_cluster(w)
    use systems, only: sysc, sys, sys_init, ok_system, remove_system
    use dynamics, only: md_relax
    use energy, only: ff_tip4p
    use utils, only: iw_text, iw_button, iw_tooltip, iw_calcwidth, iw_radiobutton,&
       iw_intstepper
    use gui_main, only: g
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use tools_io, only: string
    use param, only: kcal2ha, hartokjmol, hartoev, bohrtoa
    class(window), intent(inout), target :: w

    logical :: doquit, goodparent, ldum, hasref
    integer :: isys, nwat, tflags
    real*8 :: eb_kcal, rec_kcal, score
    type(ImVec2) :: szavail, sz0
    character(len=:,kind=c_char), allocatable, target :: str1, str2

    logical, save :: ttshown = .false. ! tooltip flag

    ! do we have a good parent window (a view)?
    goodparent = w%idparent > 0 .and. w%idparent <= nwin
    if (goodparent) goodparent = win(w%idparent)%isinit
    if (goodparent) goodparent = (win(w%idparent)%type == wintype_view)
    doquit = .not.goodparent

    ! initialize state
    if (w%firstpass) w%errmsg = ""

    if (.not.doquit) then
       ! number of water molecules
       call iw_text("Number of water molecules",highlight=.true.)
       ldum = iw_intstepper("wcnwat",w%wc_nwat,minval=int(wc_nmin,c_int),maxval=int(wc_ngen_max,c_int),&
          ndigit=3,tooltip="Number of water molecules in the cluster (2 to 100; reference energy known up to 21)")

       ! initial placement
       call igAlignTextToFramePadding()
       call iw_text("Placement",highlight=.true.)
       ldum = iw_radiobutton("Random",int=w%wc_placement,intval=0_c_int,sameline=.true.)
       call iw_tooltip("Scatter the water molecules at random positions",ttshown)
       ldum = iw_radiobutton("In a row",int=w%wc_placement,intval=1_c_int,sameline=.true.)
       call iw_tooltip("Line the water molecules up in a row",ttshown)

       ! generate a new cluster
       if (iw_button("New cluster",danger=.true.)) then
          ! stop, free, and remove the previous cluster
          if (ok_system(w%wc_isys,sys_init)) then
             sysc(w%wc_isys)%md_run = .false.
             if (sysc(w%wc_isys)%md%ready) call sysc(w%wc_isys)%md%free()
             call remove_system(w%wc_isys)
          end if
          w%wc_isys = 0
          w%errmsg = ""
          call build_water_cluster(int(w%wc_nwat),int(w%wc_placement),w%wc_isys,w%errmsg)
          w%wc_started = .false.
          w%wc_irep_text = 0
       end if
       call iw_tooltip("Generate a new cluster",ttshown)

       ! zoom buttons (touchscreens may not have a mouse wheel); only once a
       ! cluster exists and the view is displaying it
       if (ok_system(w%wc_isys,sys_init) .and. w%wc_isys == win(w%idparent)%view_selected) then
          if (iw_button("Zoom +",sameline=.true.)) call wc_zoom(0.15_c_float)
          call iw_tooltip("Zoom in",ttshown)
          if (iw_button("Zoom -",sameline=.true.)) call wc_zoom(-0.15_c_float)
          call iw_tooltip("Zoom out",ttshown)
       end if

       ! kiosk auto-start: once the generated system has finished initializing,
       ! set up TIP4P + FIRE and start the continuous relaxation
       if (ok_system(w%wc_isys,sys_init) .and. .not.w%wc_started) &
          call wc_start()

       ! scoreboard: update the on-screen text and the window status table
       isys = w%wc_isys
       if (ok_system(isys,sys_init) .and. w%wc_started) then
          if (sysc(isys)%md%ready) then
             nwat = sys(isys)%c%ncel / 3
             eb_kcal = -sysc(isys)%md%epot / kcal2ha
             hasref = (nwat >= wc_nmin .and. nwat <= wc_nmax)
             score = 0d0
             if (hasref) then
                rec_kcal = -wc_record_kjmol(nwat) / (kcal2ha * hartokjmol)
                ! score in [0,100]: 0 at zero energy (fully dissociated) and
                ! 100 at the record energy (both energies are negative)
                score = 100d0 * (sysc(isys)%md%epot*hartokjmol) / wc_record_kjmol(nwat)
                score = min(max(score,0d0),100d0)
             end if

             ! push the number + color to the on-screen scoreboard text
             call wc_update_scoreboard(isys,score,eb_kcal,hasref)

             ! window status table
             tflags = ImGuiTableFlags_None
             tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
             tflags = ior(tflags,ImGuiTableFlags_RowBg)
             tflags = ior(tflags,ImGuiTableFlags_Borders)
             tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
             str1 = "##wcstatus" // c_null_char
             sz0%x = 0._c_float
             sz0%y = 0._c_float
             if (igBeginTable(c_loc(str1),2,tflags,sz0,0._c_float)) then
                str2 = "Property" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_WidthFixed,0._c_float,0_c_int)
                str2 = "Value" // c_null_char
                call igTableSetupColumn(c_loc(str2),ImGuiTableColumnFlags_WidthFixed,0._c_float,1_c_int)
                call igTableHeadersRow()

                call status_row("Water molecules",string(nwat))
                call status_row("Binding energy (kcal/mol)",string(eb_kcal,'f',decimal=1))
                if (hasref) then
                   call status_row("Record (kcal/mol)",string(rec_kcal,'f',decimal=1))
                   call status_row("Score (0-100)",string(score,'f',decimal=1))
                end if
                if (sysc(isys)%md%nat > 0) &
                   call status_row("Max |force| (eV/A)",&
                      string(maxval(norm2(sysc(isys)%md%f,1))*hartoev/bohrtoa,'f',decimal=4))

                call igEndTable()
             end if
          end if
       end if

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

    ! quit = close the window (the run is stopped and freed in window_end)
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
    subroutine wc_zoom(ratio)
      real(c_float), intent(in) :: ratio
      call sysc(w%wc_isys)%sc%cam_zoom(ratio)
      win(w%idparent)%forcerender = .true.
    end subroutine wc_zoom

    !> Set up TIP4P + FIRE for the generated cluster and start the continuous
    !> relaxation, then add the hydrogen-bond and scoreboard representations.
    subroutine wc_start()
      use representations, only: reptype_atoms, reptype_text, reptype_axes,&
         repflavor_atoms_hbonds, repflavor_text, textpos_screen
      character(len=:), allocatable :: errmsg

      integer :: is, i

      is = w%wc_isys
      w%wc_started = .true.

      ! TIP4P relaxation engine
      sysc(is)%md_backend = ff_tip4p
      call sysc(is)%md%init(sys(is)%c,backend=ff_tip4p,mode=md_relax,errmsg=errmsg)
      if (len_trim(errmsg) > 0) then
         w%errmsg = errmsg
         sysc(is)%md_run = .false.
         return
      end if
      sysc(is)%md_run = .true.
      win(w%idparent)%forcerender = .true.

      ! hide the Cartesian axes gizmo
      do i = 1, sysc(is)%sc%nrep
         if (sysc(is)%sc%rep(i)%type == reptype_axes) &
            sysc(is)%sc%rep(i)%shown = .false.
      end do

      ! hydrogen bonds
      call sysc(is)%sc%add_representation(reptype_atoms,repflavor_atoms_hbonds)

      ! on-screen scoreboard text
      call sysc(is)%sc%add_representation(reptype_text,repflavor_text,w%wc_irep_text)
      associate (t => sysc(is)%sc%rep(w%wc_irep_text)%text%t(1))
        t%placement = textpos_screen
        t%winpos = (/0.5d0,0.90d0/)
        t%scale = 1.5d0
        t%scalewithzoom = .false.
        t%infront = .true.
      end associate
    end subroutine wc_start

    !> Write the on-screen scoreboard
    subroutine wc_update_scoreboard(is,score,eb,hasref)
      use tools_io, only: string
      integer, intent(in) :: is
      real*8, intent(in) :: score, eb
      logical, intent(in) :: hasref

      real(c_float) :: rgb(3)

      if (w%wc_irep_text < 1 .or. w%wc_irep_text > sysc(is)%sc%nrep) return
      if (sysc(is)%sc%rep(w%wc_irep_text)%text%ntext < 1) return

      ! the always-on MD run rebuilds the draw lists every frame, so the new
      ! string/color is picked up without forcing a rebuild here
      call wc_score_color(hasref,score,rgb)
      if (hasref) then
         sysc(is)%sc%rep(w%wc_irep_text)%text%t(1)%str = string(score,'f',decimal=1)
      else
         sysc(is)%sc%rep(w%wc_irep_text)%text%t(1)%str = string(eb,'f',decimal=2)
      end if
      sysc(is)%sc%rep(w%wc_irep_text)%text%t(1)%rgb = rgb
    end subroutine wc_update_scoreboard

  end subroutine draw_water_cluster

  !xx! private procedures

  !> Build a molecular system of nwat water molecules at the TIP4P reference
  !> geometry (placement: 0 = random in a box, 1 = in a row), register it as a
  !> new system, and return its id. Errors are returned in errmsg.
  subroutine build_water_cluster(nwat,placement,id,errmsg)
    use crystalseedmod, only: crystalseed
    use systems, only: add_systems_from_seeds, launch_initialization_thread
    use global, only: rborder_def
    use param, only: isformat_r_derived, bohrtoa, rad, pi
    use tools_io, only: string
    integer, intent(in) :: nwat
    integer, intent(in) :: placement
    integer, intent(inout) :: id
    character(len=:), allocatable, intent(inout) :: errmsg

    type(crystalseed), allocatable :: seed(:)
    integer, allocatable :: idlist(:)
    real*8 :: rmono(3,3), rot(3,3), roh, half, dmin, dmin2, drow, lbox, trial(3), u3(3)
    real*8, allocatable :: cen(:,:)
    integer :: i, j, k, m, iw, itry
    logical :: good, found, clash

    errmsg = ""
    if (nwat < 1) then
       errmsg = "invalid number of water molecules"
       return
    end if

    ! fresh random sequence for this cluster
    call random_seed()

    ! reference monomer (bohr): O at the origin, H's symmetric in the xy-plane
    roh = wc_roh_ang / bohrtoa
    half = 0.5d0 * wc_hoh_deg * rad
    rmono(:,1) = (/0d0, 0d0, 0d0/)
    rmono(:,2) = (/roh*cos(half),  roh*sin(half), 0d0/)
    rmono(:,3) = (/roh*cos(half), -roh*sin(half), 0d0/)

    ! monomer centers (bohr); keep O's well apart so tip4p_setup recognizes
    ! nwat distinct waters (no spurious cross-molecule bonds)
    allocate(cen(3,nwat))
    dmin = 3.0d0 / bohrtoa
    dmin2 = dmin * dmin
    if (placement == 1) then
       ! a straight row along x
       drow = 3.2d0 / bohrtoa
       do i = 1, nwat
          cen(:,i) = (/real(i-1,8)*drow, 0d0, 0d0/)
       end do
    else
       ! random positions in a cube, rejection-sampled for the minimum O-O
       ! separation. Size the box to a fixed sphere-packing fraction (~0.20)
       ! so the monomers start compact (quick clustering) but well below the
       ! random-close-packing jamming limit (~0.38) even at nwat = 100, where a
       ! tighter box would make rejection sampling fail. Grow it as a safety net.
       lbox = dmin * (real(nwat,8) * (pi/6d0) / 0.20d0)**(1d0/3d0)
       do
          good = .true.
          do i = 1, nwat
             found = .false.
             do itry = 1, 2000
                call random_number(u3)
                trial = u3 * lbox
                clash = .false.
                do j = 1, i-1
                   ! compare squared distances to avoid a sqrt per pair test
                   if (sum((trial-cen(:,j))**2) < dmin2) then
                      clash = .true.
                      exit
                   end if
                end do
                if (.not.clash) then
                   cen(:,i) = trial
                   found = .true.
                   exit
                end if
             end do
             if (.not.found) then
                good = .false.
                exit
             end if
          end do
          if (good) exit
          lbox = lbox * 1.3d0
       end do
    end if

    ! assemble the seed
    allocate(seed(1))
    seed(1)%nat = 3*nwat
    allocate(seed(1)%x(3,3*nwat),seed(1)%is(3*nwat),seed(1)%atname(3*nwat))
    seed(1)%nspc = 2
    allocate(seed(1)%spc(2))
    seed(1)%spc(1)%z = 8
    seed(1)%spc(1)%name = "O"
    seed(1)%spc(2)%z = 1
    seed(1)%spc(2)%name = "H"

    k = 0
    do iw = 1, nwat
       call wc_random_rotation(rot)
       do m = 1, 3
          k = k + 1
          seed(1)%x(:,k) = cen(:,iw) + matmul(rot,rmono(:,m))
          if (m == 1) then
             seed(1)%is(k) = 1
             seed(1)%atname(k) = "O"
          else
             seed(1)%is(k) = 2
             seed(1)%atname(k) = "H"
          end if
       end do
    end do

    ! non-periodic molecule, built in memory (no source file)
    seed(1)%ismolecule = .true.
    seed(1)%useabr = 0
    seed(1)%havesym = 0
    seed(1)%neqv = 0
    seed(1)%ncv = 0
    seed(1)%havex0 = .false.
    seed(1)%molx0 = 0d0
    seed(1)%border = rborder_def
    seed(1)%cubic = .false.
    seed(1)%isused = .true.
    seed(1)%isformat = isformat_r_derived
    seed(1)%file = ""
    seed(1)%name = "Water cluster (" // string(nwat) // ")"

    ! register the system; add_systems_from_seeds selects it in the tree, so the
    ! main view displays it
    call add_systems_from_seeds(1,seed,idlist=idlist)
    call launch_initialization_thread()
    id = 0
    if (allocated(idlist)) then
       if (size(idlist) >= 1) id = idlist(1)
    end if

  end subroutine build_water_cluster

  !> Uniform random rotation matrix. Shoemake's method builds a uniformly
  !> distributed unit quaternion; quat2mat turns it into the rotation matrix.
  subroutine wc_random_rotation(rot)
    use tools_math, only: quat2mat
    use param, only: tpi
    real*8, intent(out) :: rot(3,3)

    real*8 :: u(3), q(4)

    call random_number(u)
    q(1) = sqrt(1d0-u(1)) * sin(tpi*u(2))
    q(2) = sqrt(1d0-u(1)) * cos(tpi*u(2))
    q(3) = sqrt(u(1))     * sin(tpi*u(3))
    q(4) = sqrt(u(1))     * cos(tpi*u(3))
    rot = quat2mat(q)

  end subroutine wc_random_rotation

  !> Map the score (0 to 100) to a red-to-green color. If there is no
  !> reference energy for this cluster size (hasref = .false.) the text is
  !> plain black. Otherwise the score is stretched red -> yellow -> green.
  !> Colors are darkened so the text stays legible against the (light) scene
  !> background.
  subroutine wc_score_color(hasref,score,rgb)
    logical, intent(in) :: hasref
    real*8, intent(in) :: score
    real(c_float), intent(out) :: rgb(3)

    real(c_float), parameter :: dark = 0.72_c_float ! darken for legibility
    real*8 :: t

    ! no reference energy for this cluster size: plain black
    if (.not.hasref) then
       rgb = 0._c_float
       return
    end if

    ! score is already clamped to [0,100] by the caller
    t = score / 100d0
    if (t < 0.5d0) then
       ! red -> yellow
       rgb = (/1._c_float, real(2d0*t,c_float), 0._c_float/)
    else
       ! yellow -> green
       rgb = (/real(2d0*(1d0-t),c_float), 1._c_float, 0._c_float/)
    end if
    rgb = dark * rgb

  end subroutine wc_score_color

end submodule water_cluster
