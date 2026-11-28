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

! Routines for the atomic environments window (the ENVIRON, ECON and
! COORD keywords, for one atom).
submodule (windows) environ
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! transient markers in the view: halos around the centre, the
  ! coordinated neighbours, and the neighbours of the hovered row
  ! (color, opacity), and the radius of a halo relative to the
  ! covalent radius of the atom (the atoms are drawn at 0.7)
  real(c_float), parameter :: center_rgb(3) = (/1.0_c_float,0.6_c_float,0.1_c_float/)
  real(c_float), parameter :: coord_rgb(3) = (/0.2_c_float,0.7_c_float,1.0_c_float/)
  real(c_float), parameter :: hover_rgb(3) = (/1.0_c_float,1.0_c_float,0.2_c_float/)
  real(c_float), parameter :: halo_alpha = 0.35_c_float
  real(c_float), parameter :: hover_alpha = 0.7_c_float
  real*8, parameter :: halo_scale = 0.85d0
  real*8, parameter :: halo_radmin = 0.5d0
  ! opacity of the background of the rows with coordinated neighbours
  ! in the table of shells (in the color of their halos)
  real(c_float), parameter :: coord_bgalpha = 0.3_c_float

  ! the radial distribution: points of the curve, and height of the
  ! plot (lines)
  integer, parameter :: rdf_npts = 501
  integer, parameter :: rdf_nlines = 10
  ! margin of the neighbour search beyond what the plot needs (Å), so
  ! dragging the width does not search on every frame
  real*8, parameter :: rdf_margin = 1d0

  ! what a row of the table is (environ_state%rows; 0 = a neighbour)
  integer, parameter :: rows_shells = 1

  ! the columns of the neighbour table
  integer(c_int), parameter :: ic_shell = 0
  integer(c_int), parameter :: ic_name = 1
  integer(c_int), parameter :: ic_id = 2
  integer(c_int), parameter :: ic_cidx = 3
  integer(c_int), parameter :: ic_lvec = 4
  integer(c_int), parameter :: ic_count = 5
  integer(c_int), parameter :: ic_dist = 6
  integer(c_int), parameter :: ic_coord = 7
  integer(c_int), parameter :: ic_wecon = 8
  integer(c_int), parameter :: ic_geom = 9

contains

  !> Draw the atomic environments window: the neighbours of an atom of
  !> the system shown by the anchor view (picked in the view), in a
  !> table with a row per neighbour or per shell, as ENVIRON lists
  !> them, with the coordination (by one of three rules) and the ECoN
  !> weight of each and, per shell, the closest ideal shape to its
  !> coordinated neighbours; and the radial distribution. The atom
  !> and its neighbours are marked in the view.
  module subroutine draw_environ(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use gui_main, only: g
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_close_event,&
       iw_setpos_bottomright, iw_table_column, iw_table_headers_row, iw_calcheight,&
       iw_separatortext, iw_separatortext_height, iw_highlight_selectable, iw_cell_right,&
       iw_combo_simple, iw_intstepper
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer, parameter :: itable_neigh = 1 ! the table, for table_write_*

    logical :: doquit, goodsys, syschanged, ismol, ldum, noatoms
    integer :: isys, iview, istat, iwrite, nrow, ncol, ihover
    integer :: rj0, rj1 ! the neighbours of the radial distribution near the mouse
    real*8 :: rhover
    real(c_float) :: hbelow ! height of what is laid out below the table of neighbours
    integer :: atpick(4)
    integer(c_int) :: icol(0:9) ! position of each column in the table (-1 = not shown)
    integer(c_int) :: tflags
    real*8 :: xc(3)
    type(ImVec2) :: sz0
    type(c_ptr) :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    character(len=:,kind=c_char), allocatable, target :: str1

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)
    ! a system with no atoms has no environments
    noatoms = .false.
    if (goodsys) noatoms = (sys(isys)%c%ncel == 0)
    if (noatoms) goodsys = .false.
    ismol = .false.
    if (goodsys) ismol = sys(isys)%c%ismolecule

    ! initialize the form on the first pass
    if (w%firstpass) then
       w%env = environ_state()
       w%errmsg = ""
    end if

    ! a new system: the centre goes to the first atom
    if (goodsys) then
       if (.not.w%env%sys%is(isys) .or. syschanged) then
          call w%env%sys%set(isys)
          call cancel_pick()
          w%errmsg = ""
          call w%table_write_cancel()
          w%env%idx = (/1,0,0,0/)
          w%env%ispc = 0
          w%env%ls = environ_list()
       end if
    end if

    ! a pending pick in the view it was armed on; dropped, releasing the
    ! view, if the system is not ready or the window moved to another view
    if (w%env%picking) then
       if (goodsys .and. iview == w%env%pickview) then
          call view_pick_result(iview,w%id,isys,w%env%pick,istat,xc,atidx=atpick)
          if (istat == ipick_point .and. atpick(1) > 0) then
             w%env%idx = atpick
             call w%table_write_cancel()
          end if
          if (istat /= ipick_pending) w%env%picking = .false.
       else
          call cancel_pick()
       end if
    end if

    ! a centre that is gone after a change of geometry, or a species
    ! that is gone
    if (goodsys) then
       if (w%env%idx(1) < 1 .or. w%env%idx(1) > sys(isys)%c%ncel) w%env%idx = (/1,0,0,0/)
       if (ismol) w%env%idx(2:4) = 0
       if (w%env%ispc < 0 .or. w%env%ispc > sys(isys)%c%nspc) w%env%ispc = 0
    end if

    ihover = 0
    rhover = -1d0
    rj0 = 1
    rj1 = 0
    if (goodsys) then
       ! the neighbour list, calculated again if anything it depends on changed
       call update_list()

       ! the system
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! the centre: the atom, and its image in a crystal, in its color in the view
       call iw_text("Center",highlight=.true.,alignframe=.true.)
       ldum = draw_anchor_button(iview,isys,int(w%env%idx,c_int),"##environcenter",&
          sameline=.true.,inert=.true.)
       if (ismol) then
          call iw_tooltip("Atom whose environment is shown (marked in orange in the view): its&
             & species and its index in the molecule",ttshown)
       else
          call iw_tooltip("Atom whose environment is shown (marked in orange in the view): its&
             & species, its index in the unit cell and, if it is not in the main cell, the&
             & lattice vector of its image",ttshown)
       end if
       if (iw_button("Pick##environpick",sameline=.true.,disabled=w%env%picking)) then
          w%env%picking = .true.
          w%env%pickview = iview
          call w%env%pick%arm()
          call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the central atom",w%id)
       end if
       call iw_tooltip("Pick the atom in the view",ttshown)

       ! what is laid out below the table of neighbours, which takes
       ! the rest of the window: the radial distribution (its header
       ! and, if open, its body) and the messages
       hbelow = g%Style%ItemSpacing%y + iw_separatortext_height()
       if (w%env%secrdf%open) hbelow = hbelow + g%Style%ItemSpacing%y + rdf_body_height()
       if (sysc(isys)%md_run) hbelow = hbelow + g%Style%ItemSpacing%y + igGetTextLineHeight()
       if (len_trim(w%errmsg) > 0) hbelow = hbelow + g%Style%ItemSpacing%y + igGetTextLineHeight()

       ! the neighbours
       call iw_separatortext("Neighbors (" // string(w%env%ls%nat) // ")",w%env%secneigh)
       if (w%env%secneigh%open) call draw_neigh_section()

       ! the radial distribution
       call iw_separatortext("Radial distribution",w%env%secrdf)
       if (w%env%secrdf%open) call draw_rdf_section()

       ! the centre and its neighbours in the anchor view
       call draw_markers()

       ! the error message from writing the tables; the values wait
       ! for the end of a dynamics run
       if (sysc(isys)%md_run) call iw_text("The environment is not updated while a dynamics run&
          & is active",danger=.true.,wrap=.true.)
       if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
    end if

    if (noatoms) call iw_text("This system has no atoms",danger=.true.)

    ! bottom row: close
    call iw_setpos_bottomright(len("Close"),1)
    if (iw_button("Close")) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) then
       call cancel_pick()
       call w%end()
    end if

  contains

    !> Calculate what the window shows about the centre, each part
    !> again only if anything it depends on changed: the ECoN (the
    !> atom and the geometry), the coordination set (also the rule and
    !> the bonds), and the neighbour list (also the cutoff and the
    !> species); then which neighbours are coordinated, and the shapes
    !> of the shells.
    !> During a dynamics run nothing is calculated (the atoms move on
    !> every frame without a geometry change): the values of the atom
    !> are kept (those of another atom or system are dropped), and
    !> everything is calculated again when the run stops.
    subroutine update_list()
      use param, only: icrd_crys
      use crystalmod, only: econ_weight, econ_wmin
      real*8, allocatable :: up2d
      integer, allocatable :: up2sh, ispc0
      integer, allocatable :: iperm(:), eid(:), lvec(:,:)
      real*8, allocatable :: dist(:)
      logical, allocatable :: inec(:)
      real*8, allocatable :: econ(:), econ1(:), nd(:), nd1(:) ! the ECoN per species: not shown
      integer :: i, j, k, iat
      real*8 :: time
      logical :: newcoord, newlist

      iat = w%env%idx(1)
      associate(c => sys(isys)%c, ls => w%env%ls)
        if (sysc(isys)%md_run) then
           if (ls%ndkeysys /= w%env%sys%uid .or. ls%ndkeyiat /= iat) ls = environ_list()
           ls%mdstale = .true.
           return
        end if
        if (ls%mdstale) then
           ls%ndkeytime = -1d0
           ls%ckeytime = -1d0
           ls%keytime = -1d0
           ls%mdstale = .false.
        end if

        ! the mean weighted distance of the ECoN of the atom (all
        ! species), which gives the weights, and the neighbours with a
        ! weight of at least econ_wmin (as the integer CN of ECON); the
        ! other stages depend on the geometry too, so they follow when
        ! it changes
        time = sysc(isys)%timelastchange_geometry
        if (ls%ndkeysys /= w%env%sys%uid .or. ls%ndkeyiat /= iat .or. ls%ndkeytime /= time) then
           ls%ndkeysys = w%env%sys%uid
           ls%ndkeyiat = iat
           ls%ndkeytime = time
           allocate(econ(0:c%nspc),econ1(0:c%nspc),nd(0:c%nspc),nd1(0:c%nspc))
           call c%econ_atom(c%atcel(iat)%x,econ,econ1,nd,nd1,eid,dist,lvec)
           ls%dav = nd(0)
           allocate(inec(size(eid)))
           inec = .false.
           if (ls%dav > 0d0) inec = (econ_weight(dist,ls%dav) >= econ_wmin)
           ls%ecid = pack(eid,inec)
           ls%eclvec = lvec(:,pack((/(i,i=1,size(eid))/),inec))
        end if

        ! the coordination set
        ! (rule 1 reads the radii and bond factor of the system, which
        ! the geometry window edits before they are applied)
        time = max(sysc(isys)%timelastchange_geometry,sysc(isys)%timelastchange_rebond)
        newcoord = (ls%ckeysys /= w%env%sys%uid .or. ls%ckeyiat /= iat .or.&
           ls%ckeyrule /= w%env%cnrule .or. ls%ckeytime /= time)
        if (.not.newcoord .and. w%env%cnrule == 1) then
           newcoord = (ls%ckeybf /= sysc(isys)%bondfactor)
           if (.not.newcoord) newcoord = any(ls%ckeyrcov /= sysc(isys)%atmcov(c%spc(1:c%nspc)%z))
        end if
        if (newcoord) then
           ls%ckeysys = w%env%sys%uid
           ls%ckeyiat = iat
           ls%ckeyrule = w%env%cnrule
           ls%ckeytime = time
           ls%ckeybf = sysc(isys)%bondfactor
           ls%ckeyrcov = sysc(isys)%atmcov(c%spc(1:c%nspc)%z)
           call coordination_set(iat)
        end if

        ! the neighbour list
        newlist = .not.(ls%keysys == w%env%sys%uid .and. ls%keyiat == iat .and.&
           ls%keymode == w%env%upmode .and. ls%keyspc == w%env%ispc .and. ls%keytime == time)
        if (.not.newlist) then
           if (w%env%upmode == 0) newlist = (ls%keydist /= w%env%updist)
           if (w%env%upmode == 1) newlist = (ls%keyshell /= w%env%upshell)
        end if
        if (newlist) then
           ls%keysys = w%env%sys%uid
           ls%keyiat = iat
           ls%keymode = w%env%upmode
           ls%keyspc = w%env%ispc
           ls%keytime = time
           ls%keydist = w%env%updist
           ls%keyshell = w%env%upshell

           ! the list, around the atom in the main cell (an unallocated
           ! criterion is an absent argument)
           if (w%env%upmode == 0) then
              up2d = w%env%updist / bohrtoa
           else
              up2sh = w%env%upshell
           end if
           if (w%env%ispc > 0) ispc0 = w%env%ispc
           call c%list_near_atoms(c%atcel(iat)%x,icrd_crys,.true.,ls%nat,ls%eid,ls%dist,ls%lvec,ls%ishell,&
              up2d=up2d,up2sh=up2sh,ispc0=ispc0,nozero=.true.)

           ! within a shell (the same distance), by non-equivalent atom
           allocate(iperm(ls%nat))
           do i = 1, ls%nat
              iperm(i) = i
           end do
           do i = 2, ls%nat
              k = iperm(i)
              j = i - 1
              do while (j >= 1)
                 if (ls%ishell(iperm(j)) /= ls%ishell(k) .or.&
                    c%atcel(ls%eid(iperm(j)))%idx <= c%atcel(ls%eid(k))%idx) exit
                 iperm(j+1) = iperm(j)
                 j = j - 1
              end do
              iperm(j+1) = k
           end do
           ls%eid = ls%eid(iperm)
           ls%dist = ls%dist(iperm)
           ls%lvec = ls%lvec(:,iperm)
           ls%ishell = ls%ishell(iperm)

           ! the groups: runs in the same shell and non-equivalent atom
           if (allocated(ls%igroup)) deallocate(ls%igroup)
           allocate(ls%igroup(ls%nat+1))
           ls%ngroup = 0
           do i = 1, ls%nat
              if (i > 1) then
                 if (ls%ishell(i) == ls%ishell(i-1) .and.&
                    c%atcel(ls%eid(i))%idx == c%atcel(ls%eid(i-1))%idx) cycle
              end if
              ls%ngroup = ls%ngroup + 1
              ls%igroup(ls%ngroup) = i
           end do
           ls%igroup(ls%ngroup+1) = ls%nat + 1

           ! the ECoN weight of each neighbour
           if (allocated(ls%wecon)) deallocate(ls%wecon)
           allocate(ls%wecon(ls%nat))
           ls%wecon = 0d0
           if (ls%dav > 0d0) ls%wecon = econ_weight(ls%dist,ls%dav)
        end if

        ! the radial distribution goes with the list (calculated in
        ! its section)
        if (newlist) then
           ls%rkeysigma = -1d0
           ls%rsearch = -1d0
        end if

        ! coordinated: in the coordination set; and the shapes of the
        ! coordinated neighbours of each shell
        if (newcoord .or. newlist) then
           if (allocated(ls%iscoord)) deallocate(ls%iscoord)
           allocate(ls%iscoord(ls%nat))
           ls%iscoord = .false.
           do i = 1, ls%nat
              do j = 1, ls%ncoord
                 if (ls%cid(j) == ls%eid(i) .and. all(ls%clvec(:,j) == ls%lvec(:,i))) then
                    ls%iscoord(i) = .true.
                    exit
                 end if
              end do
           end do
           call shell_shapes()
        end if
      end associate

    end subroutine update_list

    !> The closest ideal shapes to the coordinated neighbours of each
    !> shell of the list (all its atoms in the list, whichever their
    !> kind); none for the shells without them.
    subroutine shell_shapes()
      use crystalmod, only: coordination_shape
      integer :: i, n, nsh
      logical :: last
      real*8 :: xc(3)
      real*8, allocatable :: xv(:,:)

      associate(ls => w%env%ls)
        nsh = 0
        if (ls%nat > 0) nsh = maxval(ls%ishell(1:ls%nat))
        if (allocated(ls%cshape)) deallocate(ls%cshape,ls%cshdev)
        allocate(ls%cshape(2,nsh),ls%cshdev(2,nsh),xv(3,ls%nat))
        ls%cshape = 0
        ls%cshdev = 0d0
        xc = atpos(w%env%idx(1),(/0,0,0/))

        ! the list is sorted by shell: gather the coordinated neighbours
        ! of each run
        n = 0
        do i = 1, ls%nat
           if (ls%iscoord(i)) then
              n = n + 1
              xv(:,n) = atpos(ls%eid(i),ls%lvec(:,i))
           end if
           last = (i == ls%nat)
           if (.not.last) last = (ls%ishell(i+1) /= ls%ishell(i))
           if (last) then
              if (n > 0) call coordination_shape(xc,xv(:,1:n),ls%cshape(:,ls%ishell(i)),&
                 ls%cshdev(:,ls%ishell(i)))
              n = 0
           end if
        end do
      end associate

    end subroutine shell_shapes

    !> The radial distribution around the centre (of the species of
    !> the list): the count of neighbours broadened with gaussians of
    !> unit area, and the running coordination number, up to the end
    !> of the list (beyond its last shell by four widths, in shells
    !> mode). The neighbours are searched four widths further, for the
    !> gaussians that reach into the plot, plus a margin; they are
    !> searched again only if the list changed or the plot needs them
    !> further.
    subroutine radial_distribution()
      use param, only: icrd_crys, pi
      integer, allocatable :: ispc0
      integer :: i, k, n, k0, k1, iat
      real*8 :: sig, dk, fac

      iat = w%env%idx(1)
      associate(c => sys(isys)%c, ls => w%env%ls)
        ls%rkeysigma = w%env%rdfsigma

        ! the range of the plot (Å)
        sig = w%env%rdfsigma
        if (w%env%upmode == 0) then
           ls%rmax = w%env%updist
        elseif (ls%nat > 0) then
           ls%rmax = ls%dist(ls%nat) * bohrtoa + 4d0 * sig
        else
           ls%rmax = 1d0
        end if

        ! the neighbours (Å)
        if (ls%rsearch < ls%rmax + 4d0 * sig) then
           ls%rsearch = ls%rmax + 4d0 * sig + rdf_margin
           if (w%env%ispc > 0) ispc0 = w%env%ispc
           call c%list_near_atoms(c%atcel(iat)%x,icrd_crys,.true.,ls%nrdf,ls%rid,ls%rdist,ls%rlvec,&
              up2d=ls%rsearch/bohrtoa,ispc0=ispc0,nozero=.true.)
           ls%rdist = ls%rdist * bohrtoa
        end if

        ! the broadened count, adding each gaussian where it is not
        ! negligible
        if (allocated(ls%rr)) deallocate(ls%rr,ls%rg)
        allocate(ls%rr(rdf_npts),ls%rg(rdf_npts))
        dk = ls%rmax / real(rdf_npts-1,8)
        do k = 1, rdf_npts
           ls%rr(k) = real(k-1,8) * dk
        end do
        ls%rg = 0d0
        fac = 1d0 / (sig * sqrt(2d0 * pi))
        do i = 1, ls%nrdf
           k0 = max(floor((ls%rdist(i) - 5d0*sig) / dk) + 1,1)
           k1 = min(ceiling((ls%rdist(i) + 5d0*sig) / dk) + 1,rdf_npts)
           do k = k0, k1
              ls%rg(k) = ls%rg(k) + fac * exp(-0.5d0 * ((ls%rr(k) - ls%rdist(i))/sig)**2)
           end do
        end do

        ! the running coordination number: zero at the atom, one more
        ! at each neighbour, to the end of the plot
        n = count(ls%rdist(1:ls%nrdf) <= ls%rmax)
        if (allocated(ls%sx)) deallocate(ls%sx,ls%sy)
        allocate(ls%sx(n+2),ls%sy(n+2))
        ls%sx(1) = 0d0
        ls%sy(1) = 0d0
        do i = 1, n
           ls%sx(i+1) = ls%rdist(i)
           ls%sy(i+1) = real(i,8)
        end do
        ls%sx(n+2) = ls%rmax
        ls%sy(n+2) = real(n,8)
        ls%rrefit = .true.
      end associate

    end subroutine radial_distribution

    !> Height of the body of the radial distribution section, as
    !> draw_rdf_section lays it out: the width, the plot and the line
    !> below it; the width alone if there is nothing to plot (a
    !> dynamics run that started on another atom).
    function rdf_body_height() result(h)
      real(c_float) :: h

      if (sysc(isys)%md_run .and. .not.allocated(w%env%ls%rr)) then
         h = iw_calcheight(1,0)
      else
         h = iw_calcheight(rdf_nlines+1,1)
      end if

    end function rdf_body_height

    !> The radial distribution section (the section is open): the
    !> width of the gaussians, and the plot of the broadened count and
    !> the running coordination number, calculated again if the list
    !> or the width changed (not during a dynamics run). Sets rhover
    !> to the distance under the mouse (Å), if it is on the plot, and
    !> rj0:rj1 to the neighbours within two widths of it.
    subroutine draw_rdf_section()
      character(kind=c_char,len=:), allocatable, target :: strp, strx, stry
      type(ImVec2) :: szp
      real(c_double) :: xm, ym
      integer :: n, ncn
      real*8 :: sig

      call iw_text("Width (Å)",alignframe=.true.)
      ldum = iw_dragfloat_real8("##environrdfsigma",x1=w%env%rdfsigma,speed=0.002d0,min=0.005d0,&
         max=1d0,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
      call iw_tooltip("Standard deviation of the gaussian (of unit area) that broadens each&
         & neighbor, so the area under the curve is the number of neighbors",ttshown)
      if (.not.sysc(isys)%md_run .and. w%env%ls%rkeysigma /= w%env%rdfsigma) &
         call radial_distribution()
      if (.not.allocated(w%env%ls%rr)) return

      n = size(w%env%ls%sx)
      szp%x = -1._c_float
      szp%y = iw_calcheight(rdf_nlines,0)
      strp = "##environrdfplot" // c_null_char
      if (ipBeginPlot(c_loc(strp),szp,ImPlotFlags_NoTitle)) then
         strx = "Distance (Å)" // c_null_char
         stry = "Neighbors per Å" // c_null_char
         call ipSetupAxes(c_loc(strx),c_loc(stry),ImPlotAxisFlags_None,ImPlotAxisFlags_None)
         stry = "Coordination number" // c_null_char
         call ipSetupAxis(ImAxis_Y2,c_loc(stry),ImPlotAxisFlags_AuxDefault)
         if (w%env%ls%rrefit) then
            call ipSetupAxisLimits(ImAxis_X1,0._c_double,w%env%ls%rmax,ImPlotCond_Always)
            call ipSetupAxisLimits(ImAxis_Y1,0._c_double,1.05d0*max(maxval(w%env%ls%rg),1d-10),&
               ImPlotCond_Always)
            call ipSetupAxisLimits(ImAxis_Y2,0._c_double,1.05d0*max(w%env%ls%sy(n),1d0),&
               ImPlotCond_Always)
            w%env%ls%rrefit = .false.
         end if

         call ipSetAxes(ImAxis_X1,ImAxis_Y1)
         strp = "Neighbors##environrdfcurve" // c_null_char
         call ipPlotLine(c_loc(strp),c_loc(w%env%ls%rr),c_loc(w%env%ls%rg),int(rdf_npts,c_int),&
            ImPlotLineFlags_None,0_c_int)
         call ipSetAxes(ImAxis_X1,ImAxis_Y2)
         strp = "Coordination number##environrdfcn" // c_null_char
         call ipPlotStairs(c_loc(strp),c_loc(w%env%ls%sx),c_loc(w%env%ls%sy),int(n,c_int),&
            ImPlotStairsFlags_None,0_c_int)

         ! the distance under the mouse
         if (ipIsPlotHovered()) then
            call ipGetPlotMousePos(xm,ym)
            rhover = max(real(xm,8),0d0)
         end if
         call ipEndPlot()
      end if

      ! what is under the mouse: the neighbours within two widths
      ! (rj0:rj1; the list is sorted) and the coordination number
      if (rhover >= 0d0) then
         associate(ls => w%env%ls)
           sig = w%env%rdfsigma
           do while (rj0 <= ls%nrdf)
              if (ls%rdist(rj0) >= rhover - 2d0 * sig) exit
              rj0 = rj0 + 1
           end do
           rj1 = rj0 - 1
           do while (rj1 < ls%nrdf)
              if (ls%rdist(rj1+1) > rhover + 2d0 * sig) exit
              rj1 = rj1 + 1
           end do
           ncn = count(ls%rdist(1:ls%nrdf) <= rhover)
         end associate
         call iw_text("At " // string(rhover,'f',decimal=3) // " Å: " // string(rj1-rj0+1) //&
            " neighbors within 2 widths (yellow), CN " // string(ncn))
      else
         call iw_text("Hover the plot to mark the neighbors at a distance",disabled=.true.)
      end if

    end subroutine draw_rdf_section

    !> The coordination set of cell atom iat by the rule of the form.
    !> Bonded: the bonds of the atom in the system; covalent radii:
    !> the atoms closer than the bond factor times the sum of the
    !> covalent radii (those of the system, used for its bonds); ECoN
    !> weight: the atoms with a weight of at least econ_wmin.
    subroutine coordination_set(iat)
      use param, only: icrd_crys
      integer, intent(in) :: iat
      integer :: i, n
      real*8, allocatable :: up2dsp(:,:)
      real*8 :: rcov

      associate(c => sys(isys)%c, ls => w%env%ls)
        if (w%env%cnrule == 0) then
           ! bonded
           n = 0
           if (allocated(c%nstar)) n = c%nstar(iat)%ncon
           if (allocated(ls%cid)) deallocate(ls%cid,ls%clvec)
           allocate(ls%cid(n),ls%clvec(3,n))
           do i = 1, n
              ls%cid(i) = c%nstar(iat)%idcon(i)
              ls%clvec(:,i) = c%nstar(iat)%lcon(:,i)
           end do
           ls%ncoord = n
        elseif (w%env%cnrule == 1) then
           ! covalent radii
           allocate(up2dsp(c%nspc,2))
           up2dsp(:,1) = 0d0
           rcov = sysc(isys)%atmcov(c%spc(c%atcel(iat)%is)%z)
           do i = 1, c%nspc
              up2dsp(i,2) = sysc(isys)%bondfactor * (rcov + sysc(isys)%atmcov(c%spc(i)%z))
           end do
           call c%list_near_atoms(c%atcel(iat)%x,icrd_crys,.false.,ls%ncoord,eid=ls%cid,lvec=ls%clvec,&
              up2dsp=up2dsp,nozero=.true.)
        else
           ! ECoN weight: from the ECoN stage
           ls%cid = ls%ecid
           ls%clvec = ls%eclvec
           ls%ncoord = size(ls%cid)
        end if
      end associate

    end subroutine coordination_set

    !> The neighbors section (the section is open): how far the list
    !> goes and its species, the coordination rule, and the table,
    !> with one row per neighbor or per shell.
    subroutine draw_neigh_section()
      use crystalmod, only: econ_wmin
      integer :: i
      logical :: changed
      character(len=:), allocatable :: s

      call iw_text("Up to",alignframe=.true.)
      call iw_combo_simple("##environupmode","Distance (Å)" // c_null_char // "Shells" //&
         c_null_char,w%env%upmode,sameline=.true.)
      call iw_tooltip("List the neighbors up to a distance from the atom, or up to a number&
         & of shells (atoms at the same distance)",ttshown)
      if (w%env%upmode == 0) then
         ldum = iw_dragfloat_real8("##environupdist",x1=w%env%updist,speed=0.02d0,min=0.1d0,&
            max=50d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
         call iw_tooltip("Largest distance from the atom to the neighbors in the list",ttshown)
      else
         ldum = iw_intstepper("##environupshell",w%env%upshell,minval=1_c_int,maxval=100_c_int,&
            sameline=.true.)
         call iw_tooltip("Number of shells of neighbors in the list",ttshown)
      end if

      call iw_text("Species",alignframe=.true.,gap=.true.)
      s = "All" // c_null_char
      do i = 1, sys(isys)%c%nspc
         s = s // trim(sys(isys)%c%spc(i)%name) // c_null_char
      end do
      call iw_combo_simple("##environspc",s,w%env%ispc,sameline=.true.)
      call iw_tooltip("List only the neighbors of this species",ttshown)

      call iw_text("Coordinated",alignframe=.true.)
      call iw_combo_simple("##environcnrule","Bonded" // c_null_char // "Covalent radii" //&
         c_null_char // "ECoN weight" // c_null_char,w%env%cnrule,sameline=.true.,&
         tooltips="The neighbors bonded to the atom (as the bonds of the system)" // c_null_char //&
         "The neighbors closer than the bond factor times the sum of the covalent radii&
         & (those of the system, used for its bonds)" // c_null_char //&
         "The neighbors whose weight in the ECoN of the atom is at least " //&
         string(econ_wmin,'f',decimal=1) // c_null_char,&
         ttshown=ttshown)
      call iw_tooltip("Which neighbors are coordinated to the atom: they are marked in blue in&
         & the view and, with shells as rows, in the table, which gives their geometry",ttshown)

      call iw_text("Rows",alignframe=.true.,gap=.true.)
      call iw_combo_simple("##environrows","Atoms" // c_null_char // "Shells" // c_null_char,&
         w%env%rows,sameline=.true.,changed=changed,&
         tooltips="One row per neighbor" // c_null_char //&
         "One row per shell and atom (a shell may hold several kinds of atom), with the number&
         & of atoms; the rows with coordinated neighbors are colored, and give the geometry of&
         & those in the shell" // c_null_char,&
         ttshown=ttshown)
      call iw_tooltip("What a row of the table is. The neighbors of the row under the mouse are&
         & marked in yellow in the view",ttshown)
      ! a table waiting to be written to a file is not the one shown now
      if (changed) call w%table_write_cancel()

      ! the table
      call draw_table()

    end subroutine draw_neigh_section

    !> The geometry cell of a row of the table of shells with
    !> coordinated neighbors, in shell ish: the closest ideal shape to
    !> the coordinated neighbors of the shell and its deviation, and
    !> the second closest in the tooltip. Returns the right end of the
    !> text (screen coordinates), wherever the table clips it.
    function geometry_cell(ish) result(xend)
      integer, intent(in) :: ish
      real(c_float) :: xend

      type(ImVec2) :: rmax
      character(len=:), allocatable :: str

      if (w%env%ls%cshape(1,ish) == 0) then
         ! no ideal shape with this number of neighbors
         call iw_text("-")
      else
         call iw_text(shape_label(ish,1))
         ! the text of the tooltip, only if it is going to be shown
         if (igIsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByPopup)) then
            str = "Closest ideal shape to the coordinated neighbors listed in this shell (of any&
               & kind), and the root mean square difference between their sorted angles at the&
               & atom and those of the shape (zero for a perfect match; distances do not count)"
            if (w%env%ls%cshape(2,ish) > 0) str = str // ". Next closest: " // shape_label(ish,2)
            call iw_tooltip(str,ttshown)
         end if
      end if
      call igGetItemRectMax(rmax)
      xend = rmax%x

    end function geometry_cell

    !> The closest (k = 1) or second closest (2) ideal shape to the
    !> coordinated neighbors of shell ish, and its deviation.
    function shape_label(ish,k) result(str)
      use crystalmod, only: coordination_shape_name
      integer, intent(in) :: ish, k
      character(len=:), allocatable :: str

      str = coordination_shape_name(w%env%ls%cshape(k,ish)) // " (" //&
         string(w%env%ls%cshdev(k,ish),'f',decimal=1) // "°)"

    end function shape_label

    !> The table of neighbours (the section is open): one row per
    !> neighbour, or one per group (shell and atom) if the rows are
    !> shells; then the groups with coordinated neighbours have a
    !> background and the geometry of their shell. Sets ihover to the
    !> hovered row.
    subroutine draw_table()
      use utils, only: iw_window_width
      integer :: i, i1, i2, nc
      logical :: grouped
      integer(c_int) :: coordbg
      type(ImVec2) :: pos
      real(c_float) :: xright, xover

      grouped = (w%env%rows == rows_shells)
      tflags = ImGuiTableFlags_None
      tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
      tflags = ior(tflags,ImGuiTableFlags_RowBg)
      tflags = ior(tflags,ImGuiTableFlags_Borders)
      tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
      tflags = ior(tflags,ImGuiTableFlags_ScrollY)
      ! a different table in each mode, so each keeps its own column widths
      if (grouped) then
         str1 = "##environtableshells" // c_null_char
         nrow = w%env%ls%ngroup
      else
         str1 = "##environtable" // c_null_char
         nrow = w%env%ls%nat
      end if

      ! the columns of this mode
      ncol = 0
      icol = -1
      call addcol(ic_shell)
      call addcol(ic_name)
      call addcol(ic_id)
      if (grouped) then
         call addcol(ic_count)
      else
         ! in a molecule, the cell id is the id
         if (.not.ismol) then
            call addcol(ic_cidx)
            call addcol(ic_lvec)
         end if
      end if
      call addcol(ic_dist)
      call addcol(ic_coord)
      call addcol(ic_wecon)
      if (grouped) call addcol(ic_geom)

      ! the rest of the window above what is laid out below the table,
      ! but at least a few rows (the stretch beyond them is elastic)
      call igGetCursorScreenPos(pos)
      call igGetContentRegionAvail(sz0)
      xright = pos%x + sz0%x
      xover = -huge(1._c_float)
      sz0%x = 0._c_float
      sz0%y = sz0%y - hbelow
      sz0%y = max(sz0%y,iw_calcheight(6,0,.false.))
      w%heightslack = sz0%y - iw_calcheight(6,0,.false.)
      if (.not.igBeginTable(c_loc(str1),ncol,tflags,sz0,0._c_float)) return
      ! the rows end at the vertical scroll bar, if there is one
      if (igGetScrollMaxY() > 0._c_float) xright = xright - g%Style%ScrollbarSize
      coordbg = igGetColorU32_Vec4(ImVec4(coord_rgb(1),coord_rgb(2),coord_rgb(3),coord_bgalpha))

      ! the headers, in the order of the columns
      call column(ic_shell,"Shell")
      call column(ic_name,"Atom")
      call column(ic_id,"Id")
      call column(ic_cidx,"Cell id")
      call column(ic_lvec,"Lattice vector")
      call column(ic_count,"Count")
      call column(ic_dist,"Distance (Å)")
      call column(ic_coord,"Coordinated")
      call column(ic_wecon,"ECoN weight")
      call column(ic_geom,"Geometry")
      call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)

      clipper = ImGuiListClipper_ImGuiListClipper()
      call ImGuiListClipper_Begin(clipper,nrow,igGetTextLineHeightWithSpacing())
      call w%table_write_begin(itable_neigh,iwrite,clipper)
      do while (ImGuiListClipper_Step(clipper))
         call c_f_pointer(clipper,clipper_f)
         do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
            ! the neighbours in this row: one, or the atoms of a group
            call row_range(i,i1,i2)
            nc = count(w%env%ls%iscoord(i1:i2))

            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
            if (grouped .and. nc > 0) call igTableSetBgColor(ImGuiTableBgTarget_RowBg0,coordbg,-1_c_int)
            associate(c => sys(isys)%c, ls => w%env%ls, j => w%env%ls%eid(i1))
              if (setcol(ic_shell)) &
                 call iw_cell_right(string(ls%ishell(i1)))
              ! selectable spans the row; hovering it marks its atoms in the view
              if (iw_highlight_selectable("##environrow" // string(i))) ihover = i
              if (setcol(ic_name)) &
                 call iw_text(trim(c%at(c%atcel(j)%idx)%name))
              if (setcol(ic_id)) &
                 call iw_cell_right(string(c%atcel(j)%idx))
              if (setcol(ic_count)) &
                 call iw_cell_right(string(i2-i1+1))
              if (setcol(ic_cidx)) &
                 call iw_cell_right(string(j))
              if (setcol(ic_lvec)) &
                 call iw_text(lvec_string(ls%lvec(:,i1)+w%env%idx(2:4)))
              if (setcol(ic_dist)) &
                 call iw_cell_right(string(ls%dist(i1)*bohrtoa,'f',decimal=5))
              if (setcol(ic_coord)) then
                 if (grouped) then
                    call iw_cell_right(string(nc))
                 elseif (nc > 0) then
                    call iw_text("yes")
                 end if
              end if
              if (setcol(ic_wecon)) &
                 call iw_cell_right(string(ls%wecon(i1),'f',decimal=4))
              if (nc > 0) then
                 if (setcol(ic_geom)) xover = max(xover,geometry_cell(ls%ishell(i1)) - xright)
              end if
            end associate
         end do
      end do
      call ImGuiListClipper_End(clipper)
      call ImGuiListClipper_destroy(clipper)
      call w%table_write_end(itable_neigh,table_title("Environment of atom " //&
         string(w%env%idx(1)),isys))
      call igEndTable()

      ! if the geometry (the last column) does not fit, ask the window
      ! to grow: once for each width, so the user can make it narrower
      ! afterwards
      xover = xover + g%Style%CellPadding%x + 1._c_float
      if (xover > 0._c_float) then
         xover = iw_window_width() + xover
         if (abs(xover - w%env%widthreq) >= 1._c_float) then
            w%env%widthreq = xover
            w%needwidth = xover
         end if
      end if

    end subroutine draw_table

    !> Add column id to the columns of the table (at the end).
    subroutine addcol(id)
      integer(c_int), intent(in) :: id
      icol(id) = int(ncol,c_int)
      ncol = ncol + 1
    end subroutine addcol

    !> Move to column id of the current row; false if the table does
    !> not have it, or if it is clipped.
    function setcol(id) result(ok)
      integer(c_int), intent(in) :: id
      logical :: ok
      ok = .false.
      if (icol(id) < 0) return
      ok = igTableSetColumnIndex(icol(id))
    end function setcol

    !> The header of column id, if the table has it.
    subroutine column(id,title)
      integer(c_int), intent(in) :: id
      character(len=*), intent(in) :: title
      if (icol(id) < 0) return
      call iw_table_column(title,id=id,flags=ImGuiTableColumnFlags_WidthFixed)
    end subroutine column

    !> The neighbours i1 to i2 in row irow of the table: one, or the
    !> atoms of a group if grouped.
    subroutine row_range(irow,i1,i2)
      integer, intent(in) :: irow
      integer, intent(out) :: i1, i2

      if (w%env%rows == rows_shells) then
         i1 = w%env%ls%igroup(irow)
         i2 = w%env%ls%igroup(irow+1) - 1
      else
         i1 = irow
         i2 = irow
      end if

    end subroutine row_range

    !> Mark the centre, its coordinated neighbours, and the neighbours
    !> of the hovered row in the anchor view, with translucent halos.
    subroutine draw_markers()
      use representations, only: rep_shape, shapekind_sphere, shapekind_cylinder
      type(rep_shape), allocatable :: shp(:)
      integer :: i, n, i1, i2
      real*8 :: x0(3)

      if (.not.associated(win(iview)%sc)) return
      if (win(iview)%sc%isinit == 0) return

      ! the atoms of the hovered row
      i1 = 1
      i2 = 0
      if (ihover > 0) call row_range(ihover,i1,i2)

      allocate(shp(1 + w%env%ls%ncoord + 2 * max(i2-i1+1,0) + max(rj1-rj0+1,0)))
      n = 1
      x0 = atpos(w%env%idx(1),(/0,0,0/))
      shp(1) = rep_shape(kind=shapekind_sphere,x1=x0,rad=halo_rad(w%env%idx(1)),&
         rgb=center_rgb,alpha=halo_alpha,rim=.true.)

      ! the coordinated neighbours
      do i = 1, w%env%ls%ncoord
         n = n + 1
         shp(n) = rep_shape(kind=shapekind_sphere,x1=atpos(w%env%ls%cid(i),w%env%ls%clvec(:,i)),&
            rad=halo_rad(w%env%ls%cid(i)),rgb=coord_rgb,alpha=halo_alpha,rim=.true.)
      end do

      ! the neighbours of the hovered row, with a line from the centre
      ! (it shows past the bonds)
      do i = i1, i2
         n = n + 1
         shp(n) = rep_shape(kind=shapekind_sphere,x1=atpos(w%env%ls%eid(i),w%env%ls%lvec(:,i)),&
            rad=halo_rad(w%env%ls%eid(i)),rgb=hover_rgb,alpha=hover_alpha,rim=.true.)
         n = n + 1
         shp(n) = rep_shape(kind=shapekind_cylinder,x1=x0,&
            v=reshape((/atpos(w%env%ls%eid(i),w%env%ls%lvec(:,i))-x0,0d0,0d0,0d0,0d0,0d0,0d0/),(/3,3/)),&
            rad=0.05d0,rgb=hover_rgb)
      end do

      ! the neighbours near the distance under the mouse in the plot
      if (rhover >= 0d0) then
         do i = rj0, rj1
            n = n + 1
            shp(n) = rep_shape(kind=shapekind_sphere,x1=atpos(w%env%ls%rid(i),w%env%ls%rlvec(:,i)),&
               rad=halo_rad(w%env%ls%rid(i)),rgb=hover_rgb,alpha=hover_alpha,rim=.true.)
         end do
      end if
      call win(iview)%sc%show_transient_shapes(w%uid,1,shp(1:n))

    end subroutine draw_markers

    !> Position in the scene (Cartesian, bohr; molecules: in the
    !> frame of the shapes) of cell atom iat with lattice vector lv
    !> from the centre's cell, in the image of the centre.
    function atpos(iat,lv) result(x)
      integer, intent(in) :: iat, lv(3)
      real*8 :: x(3)
      associate(c => sys(isys)%c)
        x = c%x2c(c%atcel(iat)%x + lv + w%env%idx(2:4))
        if (ismol) x = x + c%molx0
      end associate
    end function atpos

    !> Radius of the halo around cell atom iat.
    function halo_rad(iat) result(r)
      use param, only: atmcov0
      integer, intent(in) :: iat
      real*8 :: r
      associate(c => sys(isys)%c)
        r = max(halo_scale * atmcov0(c%spc(c%atcel(iat)%is)%z),halo_radmin)
      end associate
    end function halo_rad

    !> A lattice vector as text.
    function lvec_string(lv) result(str)
      integer, intent(in) :: lv(3)
      character(len=:), allocatable :: str
      str = string(lv(1)) // " " // string(lv(2)) // " " // string(lv(3))
    end function lvec_string

    !> Forget a pending pick and release the view it was armed on.
    subroutine cancel_pick()
      if (.not.w%env%picking) return
      if (w%env%pickview >= 1 .and. w%env%pickview <= nwin) &
         call win(w%env%pickview)%viewmode_release_forced(w%id)
      w%env%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_environ

end submodule environ
