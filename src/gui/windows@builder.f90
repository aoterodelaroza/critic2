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
  use types, only: vstring, neighstar, nstar_subset, substituent
  implicit none

  ! add-atoms local geometries (see addatom_template for the vertex sets)
  integer, parameter :: naddgeom = 19
  integer, parameter :: maxaddsub = 10 ! largest substituent count (pentagonal prism)
  character(len=*), parameter :: addgeom_names(naddgeom) = [character(len=24) ::&
     "atom","linear","bent","triangular","trigonal pyramid","T-shape",&
     "tetrahedral","square planar","trigonal bipyramid","square pyramid",&
     "octahedral","trigonal antiprism","pentagonal bipyramid","capped octahedron",&
     "square antiprism","trigonal dodecahedron","tricapped trigonal prism",&
     "capped square antiprism","pentagonal prism"]

  ! the fragment (substituent) library, read from flib_file on first use
  logical :: fraglib_init = .false. ! whether the library has been read
  integer :: nfraglib = 0 ! number of fragments in the library
  type(vstring), allocatable :: fraglib_name(:) ! fragment names
  integer, allocatable :: fraglib_anchor(:) ! anchor atom of each fragment
  integer, allocatable :: fraglib_attach(:) ! attachment placeholder of each fragment
  type(vstring), allocatable :: fraglib_cat(:) ! submenu each fragment belongs to
  type(vstring), allocatable :: fraglib_sub(:) ! submenu inside that one (empty = none)
  real*8, allocatable :: fraglib_radius(:) ! fictitious covalent radius, bohr (ligands; <= 0 = substituent)

  ! The most recently placed add operations, most recent first: one
  ! click on the row under the Add tools arms the same choice again.
  ! Shared by the two add tools and by every builder window, and kept
  ! for as long as the GUI runs.
  integer, parameter :: maxrecent = 10
  type recentadd
     integer :: iz = 0 ! element of the atom added (0 = this entry is a fragment)
     integer :: ig = 0 ! local geometry of the atom added
     integer :: ilib = 0 ! library index of the fragment added (0 = this entry is an atom)
  end type recentadd
  type(recentadd) :: recent(maxrecent)
  integer :: nrecent = 0 ! entries of recent(:) in use

  ! Textures of the fragment diagrams, one library fragment per slot:
  ! the entry hovered in the chooser menu (hovered one at a time), the
  ! fragment on deck in the tool panel, and one per recent entry. They
  ! are drawn in the same frame, so they cannot share a slot.
  integer, parameter :: fragslot_hover = 1
  integer, parameter :: fragslot_sel = 2
  integer, parameter :: fragslot_recent = 2 ! recent entry i uses slot fragslot_recent+i
  integer, parameter :: nfragslot = fragslot_recent + maxrecent
  integer :: fragimg_idx(nfragslot) = 0 ! library index in each slot (0 = empty)
  integer(c_int), target :: fragimg_tex(nfragslot) = 0 ! GL textures (0 = none)

  ! most neighbor directions considered when placing a ligand
  integer, parameter :: maxnbcand = 24
  real*8, parameter :: eps_dzero = 1d-10 ! degenerate-distance threshold (bohr)

  ! bond styles in the local geometry pictograms
  integer, parameter :: sty_plain = 0 ! plain line, in the plane of the diagram
  integer, parameter :: sty_wedge = 1 ! filled wedge, towards the viewer
  integer, parameter :: sty_hash = 2  ! hashed wedge, away from the viewer

contains

  ! Record an add operation at the head of the recent list: an atom of
  ! element iz with local geometry ig, or the library fragment ilib. An
  ! operation already in the list moves to the head instead of being
  ! listed twice, so the row holds the last maxrecent distinct choices.
  subroutine recent_push(iz,ig,ilib)
    integer, intent(in) :: iz, ig, ilib

    integer :: i, ipos
    type(recentadd) :: new

    if (iz <= 0 .and. ilib <= 0) return
    new%iz = iz
    new%ig = ig
    new%ilib = ilib

    ipos = 0
    do i = 1, nrecent
       if (recent(i)%iz == new%iz .and. recent(i)%ig == new%ig .and.&
          recent(i)%ilib == new%ilib) then
          ipos = i
          exit
       end if
    end do
    if (ipos == 0) then
       nrecent = min(nrecent+1,maxrecent)
       ipos = nrecent
    end if
    do i = ipos, 2, -1
       recent(i) = recent(i-1)
    end do
    recent(1) = new

  end subroutine recent_push

  !> Arm the bond-operation mode jvm on view window iview for system
  !> isys, or disarm it if it is already the active mode. Arming while
  !> another mode is active switches to the new one. Shared by the
  !> builder and the geometry window, which offer the same operations.
  module subroutine bondmode_toggle(w,iview,isys,jvm)
    use interfaces_glfw, only: glfwGetTime
    class(window), intent(inout), target :: w
    integer, intent(in) :: iview, isys, jvm

    logical :: goodview

    goodview = (iview >= 1 .and. iview <= nwin)
    if (goodview) goodview = win(iview)%isinit .and. win(iview)%isopen .and.&
       win(iview)%type == wintype_view
    if (.not.goodview) return

    if (w%builder_vm == jvm) then
       call bondmode_stop(w,iview,.true.)
    else
       if (w%builder_vm /= 0) call bondmode_stop(w,iview,.true.)
       w%builder_vm = jvm
       w%builder_isys = isys
       w%builder_time = glfwGetTime()
       ! no message: the bar shows the hint for this mode (viewmode_text)
       call win(iview)%viewmode_set_forced(jvm,idcaller=w%id)
    end if

  end subroutine bondmode_toggle

  !> Release the mode this window commanded to view window iview.
  !> goodview is whether iview is still a live view window.
  module subroutine bondmode_stop(w,iview,goodview)
    class(window), intent(inout), target :: w
    integer, intent(in) :: iview
    logical, intent(in) :: goodview

    if (goodview) call win(iview)%viewmode_release_forced(w%id,w%builder_vm)
    w%builder_vm = 0
    w%builder_isys = 0
    call w%builder_bond%clear()

  end subroutine bondmode_stop

  !> Poll the bond operation this window commanded to view window iview:
  !> release it if the view is gone, shows another system, or was taken
  !> over, and otherwise apply the edit for a delivered click. Does
  !> nothing if the active mode is not a bond operation.
  module subroutine bondmode_poll(w,iview,goodview)
    use systems, only: sys, sysc
    use interfaces_glfw, only: glfwGetTime
    class(window), intent(inout), target :: w
    integer, intent(in) :: iview
    logical, intent(in) :: goodview

    integer :: ibond(5), idxpick(4), imode, k, iord
    logical :: ok

    real*8, parameter :: stale_gap = 1d0 ! discard clicks older than this (s)

    if (.not.vm_is_bondmode(w%builder_vm)) return

    ! a geometry change since an atom was staged for a new bond
    ! invalidates its index: drop it (and its highlight)
    if (w%builder_bond%is_staged()) then
       if (w%builder_bond%is_stale(sysc(w%builder_isys)%timelastchange_geometry)) &
          call w%builder_bond%clear()
    end if

    ok = goodview
    if (ok) ok = win(iview)%viewmode == w%builder_vm .and.&
       win(iview)%vmdata%owner == w%id .and.&
       win(iview)%isys == w%builder_isys
    if (.not.ok) then
       ! the view is gone, the mode was exited (cancel key, mode combo),
       ! another window took over the pick, or the view shows another
       ! system: release the forced mode
       call bondmode_stop(w,iview,goodview)
    elseif (win(iview)%vmdata%bidx(1) > 0) then
       ! A bond click was delivered: apply the edit and stay in the mode.
       ! The stale-click guard keys on the geometry, which is right here:
       ! these edits only touch the connectivity and never renumber atoms
       ibond = win(iview)%vmdata%bidx
       win(iview)%vmdata%bidx = 0
       win(iview)%vmdata%flag = 0
       if (glfwGetTime() - w%builder_time < stale_gap .and.&
          sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) then
          ! the connectivity is the authority on whether this is a bond at
          ! all: the pick buffer can be one edit out of date, and a
          ! representation can draw contacts (van der Waals, hydrogen bonds)
          ! that were never in it
          k = sys(w%builder_isys)%c%find_bond(ibond(1),ibond(2),ibond(3:5))
          if (k == 0) then
             w%errmsg = "That contact is not a bond in the system connectivity"
          elseif (w%builder_vm == vm_builder_bondremove) then
             call sysc(w%builder_isys)%remove_bond(ibond(1),ibond(2),ibond(3:5))
          elseif (w%builder_vm == vm_builder_bondorder) then
             ! single -> double -> triple -> aromatic -> dashed -> single
             select case (sys(w%builder_isys)%c%nstar(ibond(1))%ordcon(k))
             case (1)
                iord = 2
             case (2)
                iord = 3
             case (3)
                iord = -1
             case (-1)
                iord = 0
             case default
                iord = 1
             end select
             call sysc(w%builder_isys)%set_bond_order(ibond(1),ibond(2),ibond(3:5),iord)
          end if
          if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
       end if
    elseif (win(iview)%vmdata%idx(1) > 0) then
       ! An atom click was delivered: the first click stages an atom, the
       ! second bonds the pair; clicking the staged atom again cancels it
       idxpick = win(iview)%vmdata%idx(1:4)
       imode = win(iview)%vmdata%flag
       win(iview)%vmdata%idx = 0
       if (imode == 1 .and. glfwGetTime() - w%builder_time < stale_gap .and.&
          sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) then
          if (.not.w%builder_bond%is_staged()) then
             call w%builder_bond%stage(idxpick)
          elseif (w%builder_bond%same(idxpick)) then
             call w%builder_bond%clear()
          else
             call sysc(w%builder_isys)%create_bond(w%builder_bond%idx,idxpick,&
                w%builder_vm == vm_builder_bondh,w%errmsg)
             call w%builder_bond%clear()
          end if
          ! hold the camera of the view the user is clicking in through
          ! the rebuild (the edit routines post the geometry event)
          if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
       end if
    else
       ! no click pending: stamp the reference time for the stale-click guard
       w%builder_time = glfwGetTime()
    end if

  end subroutine bondmode_poll

  !> Draw the builder window. The builder is tied to its parent view
  !> (w%idparent) and operates on the system shown in it.
  module subroutine draw_builder(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_geometry,&
       reread_system_from_file, atlisttype_ncel_frac
    use dynamics, only: md_relax
    use energy, only: ff_backend_applicable, ff_backend_default
    use gui_main, only: g, fontsize, tooltip_enabled, ColorHighlightEditDistScene
    use icons, only: icon_tex, icon_ui_editgeom, icon_ui_symmetry, icon_ui_relax
    use utils, only: iw_text, iw_button, iw_atom_button, iw_tooltip, iw_combo_simple, iw_dragfloat_real8,&
       iw_periodictable, iw_menuitem, iw_icon_togglebutton, iw_helpermark, iw_calcwidth, iw_calcheight,&
       iw_setposx_fromend, iw_close_event, iw_table_column, iw_beginmenu
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_RECALC_BONDS, BIND_NAV_MEASURE,&
       BIND_EDIT_D_A_PHI, BIND_REOPEN, BIND_PICKATOM_EXIT, BIND_PICKATOM_ALT, BIND_CANCEL
    use interfaces_glfw, only: glfwGetTime
    use tools_io, only: string, nameguess
    use tools_math, only: cross, perpendicular, axisangle2mat
    use param, only: bohrtoa, pi, eye
    class(window), intent(inout), target :: w

    integer :: isys, iview, icel, imode, nsel
    integer :: idxpick(4), nrline
    logical :: doquit, goodparent, havesys, ok, ldum, relaxing
    logical :: lstate, mdshown, lplaced
    real(c_float) :: xclick(2), xicon, hicon, yrow, reserve
    type(ImVec2) :: szchild, szrow
    character(len=:), allocatable :: errmsg, strgroup
    character(kind=c_char,len=:), allocatable, target :: strchild

    logical, save :: ttshown = .false. ! tooltip flag

    real*8, parameter :: stale_gap = 1d0 ! discard clicks older than this (s)
    ! bond order for each bond-combo index >= 1 (single, double, triple, dashed, aromatic)
    integer, parameter :: bondorder(5) = (/1,2,3,0,-1/)
    real(c_float), parameter :: palscale = 1.25_c_float ! toolbar icons, in font heights
    real(c_float), parameter :: rowspacing = 2._c_float ! vertical gap between toolbar rows (pixels)
    real(c_float), parameter :: recentscale = 2.6_c_float ! recent-entry squares, in font heights

    ! do we have a good parent window (a view)? A hidden view is still
    ! good: the builder survives the hide and the mode can be released
    iview = w%idparent
    goodparent = iview > 0 .and. iview <= nwin
    if (goodparent) goodparent = win(iview)%isinit
    if (goodparent) goodparent = (win(iview)%type == wintype_view)
    doquit = .not.goodparent

    ! initialize the state on first pass (edit_pending is NOT reset:
    ! the view sets it right before this window's first draw)
    if (w%firstpass) then
       w%errmsg = ""
       w%builder_tool = it_none
       w%builder_vm = 0
       w%builder_isys = 0
       w%builder_time = 0d0
       call w%builder_bond%clear()
       w%edit_kind = 0
       w%edit_isys = 0
       w%edit_idx = 0
       w%edit_dirty = .false.
       w%edit_time = 0d0
    end if

    ! the builder operates on the system shown in the parent view
    isys = 0
    havesys = .false.
    if (goodparent) then
       isys = win(iview)%isys
       havesys = ok_system(isys,sys_init)
    end if

    ! handle an active builder mode commanded to the parent view. The bond
    ! operations are shared with the geometry window, which offers the same
    ! three buttons, so they are polled by the shared routine
    if (vm_is_bondmode(w%builder_vm)) then
       call bondmode_poll(w,iview,goodparent)
    elseif (w%builder_vm /= 0) then
       ! a geometry change since an atom was staged for a new bond
       ! invalidates its index: drop it (and its highlight)
       if (w%builder_bond%is_staged()) then
          if (w%builder_bond%is_stale(sysc(w%builder_isys)%timelastchange_geometry)) &
             call w%builder_bond%clear()
       end if
       ok = goodparent
       if (ok) ok = win(iview)%viewmode == w%builder_vm .and.&
          win(iview)%vmdata%owner == w%id .and.&
          win(iview)%isys == w%builder_isys
       if (.not.ok) then
          ! The view is gone, the mode was exited (cancel key, mode combo),
          ! another window took over the pick, or the view shows another
          ! system: release the forced mode
          call builder_stop()
       elseif (w%builder_vm == vm_builder_addatom .or. w%builder_vm == vm_builder_addfragment) then
          ! keep the mouse ghost current
          if (w%builder_vm == vm_builder_addatom) then
             win(iview)%vmdata%tooltip_iz = w%builder_addatom_z
          elseif (allocated(w%builder_frag_name)) then
             win(iview)%vmdata%tooltip_frag = w%builder_frag_name
             win(iview)%vmdata%frag_isligand = (w%builder_frag_radius > 0d0)
          end if
          if (win(iview)%vmdata%flag /= 0) then
             ! a click was delivered: add the atoms or the fragment at the
             ! clicked position, or replace the clicked atom
             icel = win(iview)%vmdata%idx(1)
             xclick = win(iview)%vmdata%xpos
             win(iview)%vmdata%flag = 0
             win(iview)%vmdata%idx = 0
             if (glfwGetTime() - w%builder_time < stale_gap .and.&
                sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) then
                ! what was actually placed goes to the head of the recent
                ! list, so the same choice is one click away next time
                if (w%builder_vm == vm_builder_addatom) then
                   call addatom_apply(xclick,icel,lplaced)
                   if (lplaced) call recent_push(w%builder_addatom_z,w%builder_addatom_ig,0)
                else
                   call frag_place(w%builder_isys,iview,w%builder_frag_nat,&
                      w%builder_frag_z(1:w%builder_frag_nat),&
                      w%builder_frag_x(:,1:w%builder_frag_nat),w%builder_frag_ianchor,&
                      w%builder_frag_iattach,(/1d0,0d0,0d0/),w%builder_frag_radius,&
                      .true.,xclick,icel,placed=lplaced)
                   if (lplaced) call recent_push(0,0,w%builder_frag_ilib)
                end if
                if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
             end if
          else
             ! no click pending: stamp the reference time
             w%builder_time = glfwGetTime()
          end if
       elseif (win(iview)%vmdata%idx(1) > 0) then
          ! An atom click was delivered: apply the edit and stay in
          ! the mode. Discard stale clicks: those delivered while the
          ! builder was not polling.
          idxpick = win(iview)%vmdata%idx(1:4)
          imode = win(iview)%vmdata%flag
          win(iview)%vmdata%idx = 0
          if (glfwGetTime() - w%builder_time < stale_gap .and.&
             sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) then
             if (w%builder_vm == vm_builder_valence) then
                ! change valence: main pick = add a hydrogen, alternate
                ! pick = remove one
                call sysc(w%builder_isys)%change_valence(idxpick(1),imode)
             elseif (w%builder_vm == vm_builder_remove .and. imode == 1) then
                ! remove atoms (main pick only): the atom and its
                ! terminal hydrogens
                call sysc(w%builder_isys)%remove_atom_hydrogens(idxpick(1))
             elseif (w%builder_vm == vm_builder_trim .and. imode == 1) then
                ! trim branch (main pick only): the same, plus every
                ! piece the removal disconnects but the main one
                call sysc(w%builder_isys)%trim_branch(idxpick(1))
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

    ! number of measure-selected atoms in the parent view: the geometry
    ! editor works on them
    nsel = 0
    if (goodparent) then
       if (associated(win(iview)%sc)) nsel = win(iview)%sc%nmsel
    end if

    ! the keybinding request toggles an edit session: two selected
    ! atoms toggle edit distance, three edit angle, four edit dihedral,
    ! and otherwise deactivate the active session, if any. Either way
    ! the geometry editor becomes the tool on display, so the panel
    ! shows the session this opened.
    if (w%edit_pending .or. (w%focused() .and. is_bind_event(BIND_EDIT_D_A_PHI))) then
       w%edit_pending = .false.
       if (nsel >= 2 .and. nsel <= 4) then
          call edit_toggle(nsel)
          w%builder_tool = it_edit
       elseif (w%edit_kind /= 0) then
          call w%edit_stop()
          w%builder_tool = it_edit
       end if
    end if

    ! per-frame guards on the active edit session
    if (w%edit_kind /= 0) then
       if (edit_session_ok()) then
          w%edit_time = glfwGetTime()
       else
          call w%edit_stop()
       end if
    end if

    ! No tool is on display unless one is running
    if (w%builder_vm == 0 .and. w%edit_kind == 0) w%builder_tool = it_none

    ! header: the system this builder operates on
    call iw_text("System",highlight=.true.)
    if (havesys) then
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)
    else
       call iw_text("none",disabled=.true.,sameline=.true.)
    end if

    ! the tools, one labelled row per group
    xicon = iw_calcwidth(10,0)
    hicon = palscale * fontsize%y + 2._c_float * g%Style%FramePadding%y
    szrow%x = g%Style%ItemSpacing%x
    szrow%y = rowspacing
    call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,szrow)

    call row_label("Add")
    call tool_button(vm_builder_addatom,.true.)
    call tool_button(vm_builder_addfragment,.false.)
    call recent_row()

    call row_label("Valence")
    call tool_button(vm_builder_valence,.true.)

    call row_label("Remove")
    call tool_button(vm_builder_remove,.true.)
    call tool_button(vm_builder_trim,.false.)

    call row_label("Bonds")
    call tool_button(vm_builder_bond,.true.)
    call tool_button(vm_builder_bondh,.false.)
    call tool_button(vm_builder_bondremove,.false.)
    call tool_button(vm_builder_bondorder,.false.)

    ! the geometry editor is driven by the atoms selected in the view, so
    ! this row carries the instructions
    call row_label("Edit",help="Select two, three, or four atoms in the view ("//&
       trim(get_bind_keyname(BIND_NAV_MEASURE))//"), then press this button or "//&
       trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//": two atoms edit a distance, three an"//&
       " angle, and four a dihedral.")
    call tool_button(it_edit,.true.)

    call row_label("Symmetry")
    call row_cursor(.true.)
    if (iw_icon_togglebutton("##buildersymm",icon_tex(icon_ui_symmetry),"Sy",&
       disabled=.not.havesys,scale=palscale)) then
       w%errmsg = ""
       call sysc(isys)%refine_symmetry(w%errmsg)
       ! hold the camera on the parent view's scene (the main scene when
       ! the parent is the main view, its private scene otherwise)
       if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
    end if
    call iw_tooltip("Refine symmetry: move the atoms onto their symmetry positions exactly",&
       ttshown,whendisabled=.true.)

    ! the current group, after the tool: space group in a crystal, point
    ! group in a molecule
    if (havesys) then
       strgroup = "n/a"
       if (sys(isys)%c%ismolecule) then
          if (sys(isys)%c%pg%avail) strgroup = trim(sys(isys)%c%pg%symbol)
       else
          if (sys(isys)%c%spgavail) strgroup = trim(sys(isys)%c%spg%international_symbol)
       end if
       call igSameLine(0._c_float,-1._c_float)
       call igSetCursorPosY(yrow + 0.5_c_float * (hicon - igGetTextLineHeight()))
       call iw_text("(" // strgroup // ")")
    end if

    ! the relaxation: the button starts and stops the run, and the
    ! method and the convergence threshold follow it on the same row
    relaxing = .false.
    mdshown = .false.
    if (havesys) then
       relaxing = sysc(isys)%md_run .and. sysc(isys)%md%mode == md_relax
       mdshown = sysc(isys)%md%ready .and. sysc(isys)%md%nat > 0 .and.&
          sysc(isys)%md%mode == md_relax
       ! resolve the method before the button (the combo below would
       ! only snap an inapplicable backend after this frame's click)
       if (sysc(isys)%md_backend < 0 .or.&
          .not.ff_backend_applicable(sysc(isys)%md_backend,sys(isys)%c)) &
          sysc(isys)%md_backend = ff_backend_default(sys(isys)%c)
    end if
    call row_label("Relax")
    call row_cursor(.true.)
    lstate = relaxing
    ok = havesys
    if (ok) ok = relaxing .or. .not.sysc(isys)%md_run
    if (iw_icon_togglebutton("##builderrelax",icon_tex(icon_ui_relax),"Rx",state=lstate,&
       disabled=.not.ok,danger=relaxing,scale=palscale)) then
       if (relaxing) then
          call sysc(isys)%md_stop()
          if (associated(win(iview)%sc)) win(iview)%forcerender = .true.
       else
          ! release any edit session (a running relaxation ends it anyway)
          call w%edit_stop()
          call sysc(isys)%md_start(md_relax,errmsg)
          w%errmsg = errmsg
          if (len_trim(w%errmsg) == 0) then
             if (associated(win(iview)%sc)) win(iview)%forcerender = .true.
          end if
       end if
    end if
    if (relaxing) then
       call iw_tooltip("Stop the geometry relaxation ("//&
          trim(get_bind_keyname(BIND_CANCEL))//")",ttshown)
    else
       call iw_tooltip("Relax the geometry to the nearest energy minimum, stopping when"//&
          " the maximum force falls below the threshold",ttshown,whendisabled=.true.)
    end if
    if (havesys) then
       ! method combo and force convergence threshold, next to the button
       call draw_ff_backend_combo(isys,"##builderrelaxmethod",15)
       call iw_tooltip("Method for the calculation of energies and forces",ttshown)
       ldum = iw_dragfloat_real8("Fmax (eV/Å)##relaxfconv",x1=sysc(isys)%md%fconv,&
          speed=0.0005d0,min=0.0001d0,max=1d0,decimal=4,sameline=.true.,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Stop the relaxation when the maximum force falls below this value",ttshown)

       ! live status: energy and maximum force
       if (mdshown) then
          call iw_text("Energy: "//string(sysc(isys)%md%epot,'f',decimal=6)//" Ha")
          call iw_text("Max force: "//string(sysc(isys)%md%maxforce(),'f',decimal=4)//" eV/Å")
          if (.not.sysc(isys)%md_run .and. sysc(isys)%md%converged()) &
             call iw_text("(converged)",highlight=.true.,sameline=.true.)
       end if

       ! surface force-evaluation errors from the run
       if (allocated(sysc(isys)%md%errmsg)) then
          if (len_trim(sysc(isys)%md%errmsg) > 0) &
             call iw_text(trim(sysc(isys)%md%errmsg),danger=.true.)
       end if
    end if
    call igPopStyleVar(1)

    ! one text line of air between the tool rows and the buttons below
    szrow%x = 0._c_float
    szrow%y = fontsize%y
    call igDummy(szrow)

    ! the two actions that discard work: the bonds the user made by hand,
    ! and every edit in the session
    if (iw_button("Rebond",danger=.true.,disabled=.not.havesys)) then
       w%errmsg = ""
       call sysc(isys)%rebond()
    end if
    call iw_tooltip("Recompute the bond connectivity for this system, discarding the bonds"//&
       " created or removed by hand ("//trim(get_bind_keyname(BIND_RECALC_BONDS))//")",&
       ttshown,whendisabled=.true.)
    if (iw_button("Restore",danger=.true.,disabled=.not.havesys,sameline=.true.)) then
       w%errmsg = ""
       ! drop any edit session without rebuilding (the geometry is
       ! discarded) and treat the system as not ready for the rest of this
       ! frame: the re-read deallocates sys(isys)%c until the init thread runs
       w%edit_dirty = .false.
       call w%edit_stop()
       call reread_system_from_file(isys)
       havesys = .false.
    end if
    call iw_tooltip("Restore the system to the original geometry it had when it was"//&
       " first opened, discarding every edit ("//&
       trim(get_bind_keyname(BIND_REOPEN))//")",ttshown,whendisabled=.true.)

    ! shortcuts to the windows where the structure can be inspected in
    ! detail: the tabs of View/Edit Geometry, and Dynamics
    call window_button("Atoms",geomtab_atoms)
    call window_button("Cell",geomtab_cell)
    call window_button("Bonds",geomtab_bonds)
    call window_button("Symmetry",geomtab_symmetry)
    call window_button("Dynamics",geomtab_none)

    ! the options of the tool on display, in a child region so that only
    ! this part scrolls: the footer buttons stay in view
    call igSeparator()
    nrline = 0
    if (len_trim(w%errmsg) > 0) nrline = 1
    reserve = iw_calcheight(1,nrline,.true.) + 2._c_float*g%Style%ItemSpacing%y
    strchild = "##buildertoolpanel" // c_null_char
    szchild%x = 0._c_float
    szchild%y = -reserve
    if (igBeginChild_Str(c_loc(strchild),szchild,.false._c_bool,ImGuiWindowFlags_None)) then
       select case (w%builder_tool)
       case (vm_builder_addatom)
          call panel_pick(vm_builder_addatom)
          call draw_addatom_geom_grid(w%builder_addatom_ig,1.9_c_float)
       case (it_edit)
          call panel_edit()
       case (it_none)
          call panel_text("Choose a tool from the panel above to edit the structure. Press "//&
             trim(get_bind_keyname(BIND_CANCEL))//" to cancel the edit.")
       case default
          call panel_pick(w%builder_tool)
       end select
    end if
    call igEndChild()
    call igSeparator()

    ! transient highlight of the latched atoms
    if (w%edit_kind /= 0) &
       call sysc(w%edit_isys)%highlight_atoms(.true.,w%edit_idx(1,1:w%edit_kind),&
       atlisttype_ncel_frac,spread(ColorHighlightEditDistScene,2,w%edit_kind))

    ! transient highlight of the atom staged for a new bond
    if (w%builder_bond%is_staged()) &
       call sysc(w%builder_isys)%highlight_atoms(.true.,w%builder_bond%idx(1:1),&
       atlisttype_ncel_frac,spread(ColorHighlightEditDistScene,2,1))

    ! error message, if any
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! close button and binds
    call iw_setposx_fromend(5,1)
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit the window (window_end releases the forced mode if still active)
    if (doquit) &
       call w%end()

  contains
    ! Open a row of the toolbar with its section label, latching the top
    ! of the row in yrow. The label is centered on the icon buttons,
    ! which are taller than a text line. help adds a help marker after
    ! the label, for a row that needs instructions.
    ! A button that opens the window where the structure can be inspected
    ! in detail: tab itab of the View/Edit Geometry window, or (for
    ! geomtab_none) the Dynamics window. The cell tab exists only in
    ! crystals.
    subroutine window_button(str,itab)
      character(len=*), intent(in) :: str
      integer, intent(in) :: itab

      integer :: idum
      logical :: ok

      ok = havesys
      if (ok .and. itab == geomtab_cell) ok = .not.sys(isys)%c%ismolecule
      if (iw_button(str//"##builderwin"//str,disabled=.not.ok,sameline=.true.)) then
         if (itab == geomtab_none) then
            idum = stack_create_window(wintype_dynamics,.true.,idparent=iview,orraise=-1)
         else
            idum = stack_create_window(wintype_geometry,.true.,isys=isys,orraise=-1)
            if (idum > 0) win(idum)%geometry_seltab = itab
         end if
      end if
      if (itab == geomtab_none) then
         call iw_tooltip("Open the Dynamics window: molecular dynamics and interactive"//&
            " manipulation of this system",ttshown,whendisabled=.true.)
      else
         call iw_tooltip("Open the "//str//" tab of the View/Edit Geometry window for this"//&
            " system",ttshown,whendisabled=.true.)
      end if

    end subroutine window_button

    subroutine row_label(str,help)
      character(len=*), intent(in) :: str
      character(len=*), intent(in), optional :: help

      yrow = igGetCursorPosY()
      call igSetCursorPosY(yrow + 0.5_c_float * (hicon - igGetTextLineHeight()))
      call iw_text(str,highlight=.true.)
      if (present(help)) call iw_helpermark(help)

    end subroutine row_label

    ! Place the cursor for the next widget of the current toolbar row: at
    ! the icon column for the first one, after the previous one
    ! otherwise, and always at the top of the row. Without the last part
    ! they would line up with the label, which sits lower to center it on
    ! the buttons (igSameLine restores the y of the previous line, not of
    ! the previous item).
    subroutine row_cursor(first)
      logical, intent(in) :: first

      if (first) then
         call igSameLine(xicon,0._c_float)
      else
         call igSameLine(0._c_float,-1._c_float)
      end if
      call igSetCursorPosY(yrow)

    end subroutine row_cursor

    ! Draw the toolbar button of tool itool, highlighted while the tool
    ! is armed on the view. first is whether it opens its row.
    subroutine tool_button(itool,first)
      integer, intent(in) :: itool
      logical, intent(in) :: first

      logical :: state, disabled, lpop, ldum_
      integer :: iz
      integer(c_int) :: itex
      character(len=:), allocatable :: fallback

      call row_cursor(first)
      if (itool == it_edit) then
         ! the geometry editor needs the atoms selected in the view
         state = (w%edit_kind /= 0)
         disabled = .not.(w%edit_kind /= 0 .or. (havesys .and. nsel >= 2 .and. nsel <= 4))
         itex = icon_tex(icon_ui_editgeom)
         fallback = "dA"
      else
         state = (w%builder_vm == itool)
         disabled = .not.havesys
         call viewmode_icon(itool,itex,fallback)
      end if
      if (itool == vm_builder_addatom .or. itool == vm_builder_addfragment) then
         ! the button opens the list of what the tool adds, the periodic
         ! table or the fragment library, whether or not it is armed
         ldum_ = iw_icon_togglebutton("##bldtool"//string(itool),itex,fallback,&
            state=state,disabled=disabled,scale=palscale,popupcontext=lpop,&
            popupflags=ImGuiPopupFlags_MouseButtonLeft)
         if (lpop) then
            if (itool == vm_builder_addatom) then
               iz = iw_periodictable()
               if (iz > 0) then
                  w%errmsg = ""
                  w%builder_addatom_z = iz
                  w%builder_addatom_ig = addatom_prefgeom(iz)
                  w%builder_tool = vm_builder_addatom
                  if (w%builder_vm /= vm_builder_addatom) &
                     call builder_toggle(vm_builder_addatom)
                  call igCloseCurrentPopup()
               end if
            else
               ! the menu entries arm the tool themselves, on the way out
               call fragment_library_menu()
            end if
            call igEndPopup()
         end if
      elseif (iw_icon_togglebutton("##bldtool"//string(itool),itex,fallback,&
         state=state,disabled=disabled,scale=palscale)) then
         call pal_click(itool)
      end if
      ! the tooltip text is only built when it is about to be shown: a
      ! palette of these would rebuild ten of them every frame
      if (tooltip_enabled) then
         if (igIsItemHovered(ior(ImGuiHoveredFlags_AllowWhenBlockedByPopup,&
            ImGuiHoveredFlags_AllowWhenDisabled))) &
            call iw_tooltip(tool_tooltip(itool),ttshown,whendisabled=.true.)
      end if

    end subroutine tool_button

    ! The row of recently placed atoms and fragments, under the Add
    ! tools: each entry redraws the choice it stands for (the local
    ! geometry with the element symbol, or the fragment diagram) and one
    ! click arms the tool with it again. Nothing is drawn until
    ! something has been placed.
    subroutine recent_row()
      integer :: i
      real(c_float) :: side, xrun, xmax
      type(ImVec2) :: sz, p0, p1, szavail
      type(c_ptr) :: dl
      integer(c_int) :: col
      logical :: hovered
      character(len=:,kind=c_char), allocatable, target :: str1

      if (nrecent == 0) return
      side = recentscale * fontsize%y
      sz%x = side
      sz%y = side
      dl = igGetWindowDrawList()
      str1 = "##recentadd" // c_null_char
      call igGetContentRegionAvail(szavail)
      xmax = igGetCursorPosX() + szavail%x
      xrun = xicon
      do i = 1, nrecent
         ! first entry, or one that no longer fits: start a row at the
         ! icon column; the previous entry has already ended the line
         if (i == 1 .or. xrun + g%Style%ItemSpacing%x + side > xmax) then
            call igSetCursorPosX(xicon)
            yrow = igGetCursorPosY()
            xrun = xicon
         else
            call row_cursor(.false.)
            xrun = xrun + g%Style%ItemSpacing%x
         end if
         xrun = xrun + side

         call igPushID_Int(int(i,c_int))
         if (igInvisibleButton(c_loc(str1),sz,ImGuiButtonFlags_None)) call recent_click(i)
         hovered = igIsItemHovered(ImGuiHoveredFlags_None)
         call igGetItemRectMin(p0)
         p1%x = p0%x + side
         p1%y = p0%y + side
         if (hovered) then
            col = igGetColorU32_Col(ImGuiCol_ButtonHovered,0.6_c_float)
            call ImDrawList_AddRectFilled(dl,p0,p1,col,0._c_float,0_c_int)
         end if
         col = igGetColorU32_Col(ImGuiCol_Border,1._c_float)
         call ImDrawList_AddRect(dl,p0,p1,col,0._c_float,0_c_int,1._c_float)
         call recent_paint(i,p0,p1,side,dl)
         if (hovered) call iw_tooltip(recent_label(i))
         call igPopID()
      end do

    end subroutine recent_row

    ! Draw what recent entry i stands for inside the square p0-p1.
    subroutine recent_paint(i,p0,p1,side,dl)
      integer, intent(in) :: i
      type(ImVec2), intent(in) :: p0, p1
      real(c_float), intent(in) :: side
      type(c_ptr), intent(in) :: dl

      integer :: islot
      integer(c_int) :: itex, col
      type(ImVec2) :: uv0, uv1
      character(len=:), allocatable :: fallback

      if (recent(i)%ilib > 0) then
         ! a fragment: its diagram, or the generic glyph if it has none
         islot = fragslot_recent + i
         if (fragimg_ensure(recent(i)%ilib,islot)) then
            uv0%x = 0._c_float
            uv0%y = 0._c_float
            uv1%x = 1._c_float
            uv1%y = 1._c_float
            col = igGetColorU32_Col(ImGuiCol_Text,1._c_float)
            call ImDrawList_AddImage(dl,int(fragimg_tex(islot),c_intptr_t),p0,p1,uv0,uv1,col)
         else
            call viewmode_icon(vm_builder_addfragment,itex,fallback)
            if (itex /= 0) then
               uv0%x = 0._c_float
               uv0%y = 0._c_float
               uv1%x = 1._c_float
               uv1%y = 1._c_float
               col = igGetColorU32_Col(ImGuiCol_Text,1._c_float)
               call ImDrawList_AddImage(dl,int(itex,c_intptr_t),p0,p1,uv0,uv1,col)
            end if
         end if
      else
         ! an atom: the local geometry, with the element at its center
         call addatom_geom_paint(recent(i)%ig,p0,side,dl,iz=recent(i)%iz)
      end if

    end subroutine recent_paint

    ! The tooltip of recent entry i: the choice it arms.
    function recent_label(i) result(str)
      integer, intent(in) :: i
      character(len=:), allocatable :: str

      str = ""
      if (recent(i)%ilib > 0) then
         if (recent(i)%ilib <= nfraglib) &
            str = "Add " // trim(fraglib_name(recent(i)%ilib)%s) // " again"
      else
         str = "Add " // trim(nameguess(recent(i)%iz,.true.)) // " (" //&
            trim(addgeom_names(recent(i)%ig+1)) // ") again"
      end if

    end function recent_label

    ! Arm the tool that recent entry i stands for, with the same element
    ! and local geometry, or the same fragment. The tool already armed
    ! only changes its choice: this is not a toggle.
    subroutine recent_click(i)
      integer, intent(in) :: i

      w%errmsg = ""
      if (recent(i)%ilib > 0) then
         if (.not.addfrag_load(recent(i)%ilib)) return
         w%builder_tool = vm_builder_addfragment
         if (w%builder_vm /= vm_builder_addfragment) &
            call builder_toggle(vm_builder_addfragment)
      else
         w%builder_addatom_z = recent(i)%iz
         w%builder_addatom_ig = recent(i)%ig
         w%builder_tool = vm_builder_addatom
         if (w%builder_vm /= vm_builder_addatom) &
            call builder_toggle(vm_builder_addatom)
      end if

    end subroutine recent_click

    ! Handle a click on the palette button of tool itool: show its
    ! options in the panel, and arm it (or disarm it, if it was the
    ! active tool).
    subroutine pal_click(itool)
      integer, intent(in) :: itool

      w%errmsg = ""
      w%builder_tool = itool
      if (itool == it_edit) then
         ! the same dispatch as the keybinding, so the button and the key
         ! never disagree: the selection chooses the kind, and asking for
         ! the kind already running stops it
         if (nsel >= 2 .and. nsel <= 4) then
            call edit_toggle(nsel)
         elseif (w%edit_kind /= 0) then
            call w%edit_stop()
         end if
      else
         ! the two add tools never reach this: their buttons open their
         ! chooser, and the pick there arms them
         call builder_toggle(itool)
      end if

    end subroutine pal_click

    ! The tooltip for the palette button of tool itool: its name and
    ! what clicking in the view will do.
    function tool_tooltip(itool) result(str)
      integer, intent(in) :: itool
      character(len=:), allocatable :: str

      character(len=:), allocatable :: descr

      if (itool == it_edit) then
         str = "Edit geometry: change the distance, angle, or dihedral of the atoms"//&
            " selected in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//"):"//&
            " two atoms for a distance, three for an angle, four for a dihedral (shortcut: "//&
            trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//")"
      else
         str = trim(vmnames(itool))
         descr = mode_descr(itool)
         if (len_trim(descr) > 0) str = str // ": " // trim(descr)
      end if

    end function tool_tooltip

    ! What clicking in the view does in mode imode, as the parent view
    ! itself would announce it (empty if there is no parent view).
    function mode_descr(jmode) result(descr)
      integer, intent(in) :: jmode
      character(len=:), allocatable :: descr

      character(len=:), allocatable :: hint, picklbl, altlbl

      descr = ""
      if (goodparent) &
         call viewmode_text(win(iview),jmode,hint,descr,picklbl,altlbl)

    end function mode_descr

    ! A paragraph of text, wrapped to the width of the panel.
    subroutine panel_text(str,faded)
      character(len=*), intent(in) :: str
      logical, intent(in), optional :: faded

      logical :: faded_
      character(len=:,kind=c_char), allocatable, target :: strl

      faded_ = .false.
      if (present(faded)) faded_ = faded
      strl = trim(str) // c_null_char
      if (faded_) &
         call igPushStyleColor_Vec4(ImGuiCol_Text,g%Style%Colors(ImGuiCol_TextDisabled+1))
      call igPushTextWrapPos(0._c_float)
      call igTextWrapped(c_loc(strl))
      call igPopTextWrapPos()
      if (faded_) call igPopStyleColor(1)

    end subroutine panel_text

    ! The ways out of mode itool
    function mode_exits(itool) result(str)
      integer, intent(in) :: itool
      character(len=:), allocatable :: str

      str = "Stop with " // trim(get_bind_keyname(BIND_CANCEL))
      if (itool == vm_builder_valence) then
         str = str // ", with " // trim(get_bind_keyname(BIND_PICKATOM_ALT)) //&
            " away from an atom"
      else
         str = str // ", with " // trim(get_bind_keyname(BIND_PICKATOM_ALT))
      end if
      if (vm_exits_on_empty(itool)) &
         str = str // ", with " // trim(get_bind_keyname(BIND_PICKATOM_EXIT)) //&
         " on empty space"
      if (itool == vm_builder_addatom .or. itool == vm_builder_addfragment) then
         ! their buttons reopen the chooser instead of disarming
         str = str // ", or by choosing another tool."
      else
         str = str // ", by clicking the tool button again, or by choosing another tool."
      end if

    end function mode_exits

    ! Panel for a tool whose only interaction is clicking in the view:
    ! what the click does, and how far along the two-click bond
    ! operations are.
    subroutine panel_pick(itool)
      integer, intent(in) :: itool

      integer :: iat
      logical :: okat
      character(len=:), allocatable :: descr

      descr = mode_descr(itool)
      if (len_trim(descr) > 0) call panel_text(descr)
      call panel_text(mode_exits(itool),faded=.true.)
      if (vm_is_bondmode(itool) .and. w%builder_bond%is_staged()) then
         iat = w%builder_bond%idx(1)
         okat = havesys
         if (okat) okat = (iat >= 1 .and. iat <= sys(isys)%c%ncel)
         if (okat) &
            call iw_text("First atom: "//&
            trim(sys(isys)%c%at(sys(isys)%c%atcel(iat)%idx)%name)//string(iat)//&
            " (click the second one)",highlight=.true.)
      end if

    end subroutine panel_pick

    ! The fragment library as nested menus: one submenu per category, in
    ! order of first appearance, with the fragments of a "parent/child"
    ! category one level deeper.
    subroutine fragment_library_menu()

      integer :: i, j, k

      call fraglib_ensure()
      if (nfraglib == 0) then
         call iw_text("No fragments available")
         return
      end if
      do k = 1, nfraglib
         ! skip a category whose submenu has already been drawn
         if (.not.fraglib_isfirst(k,.false.)) cycle
         if (iw_beginmenu(trim(fraglib_cat(k)%s))) then
            ! the fragments listed directly under this category
            do i = 1, nfraglib
               if (fraglib_cat(i)%s /= fraglib_cat(k)%s) cycle
               if (len_trim(fraglib_sub(i)%s) > 0) cycle
               call fragment_menuitem(i)
            end do
            ! then one submenu per subcategory
            do j = 1, nfraglib
               if (fraglib_cat(j)%s /= fraglib_cat(k)%s) cycle
               if (len_trim(fraglib_sub(j)%s) == 0) cycle
               if (.not.fraglib_isfirst(j,.true.)) cycle
               if (iw_beginmenu(trim(fraglib_sub(j)%s))) then
                  do i = 1, nfraglib
                     if (fraglib_cat(i)%s /= fraglib_cat(k)%s) cycle
                     if (fraglib_sub(i)%s /= fraglib_sub(j)%s) cycle
                     call fragment_menuitem(i)
                  end do
                  call igEndMenu()
               end if
            end do
            call igEndMenu()
         end if
      end do

    end subroutine fragment_library_menu

    ! Panel for the geometry editor: with no session running, how to
    ! start one; with a session running, the controls for the distance,
    ! angle, or dihedral being edited (the kind follows from how many
    ! atoms were selected).
    subroutine panel_edit()

      integer :: iside, ibold, ibnew
      integer :: nmode, imode, islot, imodes(4)
      real*8 :: dist, ang
      character(len=:), allocatable :: stropt

      if (w%edit_kind == 0) then
         call panel_text("Select two atoms in the view ("//&
            trim(get_bind_keyname(BIND_NAV_MEASURE))//") to edit a distance, three for an"//&
            " angle, or four for a dihedral, then press "//&
            trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//" or the tool button.")
         if (nsel > 0) &
            call iw_text(string(nsel)//" atoms selected in the view",highlight=.true.)
         return
      end if

      ! the latched atoms, common to the three kinds
      call edit_atom_labels()

      select case (w%edit_kind)
      case (2)
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
            if (w%edit_fragok(iside)) stropt = stropt // "Translate fragment"//c_null_char
            call iw_combo_simple("Atom "//string(iside)//"##editdistmove"//string(iside),&
               stropt,w%edit_imove(iside))
            call iw_tooltip("What moves when the distance changes: nothing (fixed), the atom,"//&
               " or the whole fragment attached to it",ttshown)
         end do

         ! distance drag-float (bohr internally, shown in angstrom, no upper bound)
         dist = norm2(edit_pos(2) - edit_pos(1))
         if (edit_drag("Distance (Å)##editdistdrag",dist,0.01d0,3,0.1d0,scal=bohrtoa)) &
            call editdist_apply(dist)
         call iw_tooltip("Distance between the two atoms (Å)",ttshown)
         call edit_apply_button("","distance")
      case (3)
         ! central atom combo; while the central atom moves, the
         ! terminals are held fixed (rotation needs a fixed vertex)
         stropt = "Fixed"//c_null_char//"Translate atom"//c_null_char
         if (w%edit_fragok(2)) stropt = stropt // "Translate group"//c_null_char
         call iw_combo_simple("Atom 2 (center)##editangmove2",stropt,w%edit_imove(2))
         call iw_tooltip("What moves when the angle changes: nothing (fixed), the central"//&
            " atom, or the whole group attached to it",ttshown)
         if (w%edit_imove(2) > 0) then
            w%edit_imove(1) = 0
            w%edit_imove(3) = 0
         end if

         ! terminal atom combos (rotation available only with a fixed center)
         do iside = 1, 3, 2
            stropt = "Fixed"//c_null_char//"Rotate atom"//c_null_char
            if (w%edit_fragok(iside)) stropt = stropt // "Rotate group"//c_null_char
            if (w%edit_imove(2) > 0) call igBeginDisabled(.true._c_bool)
            call iw_combo_simple("Atom "//string(iside)//"##editangmove"//string(iside),&
               stropt,w%edit_imove(iside))
            if (w%edit_imove(2) > 0) call igEndDisabled()
            call iw_tooltip("What moves when the angle changes: nothing (fixed), the"//&
               " terminal atom, or the whole group attached to it, rotating about the"//&
               " central atom",ttshown)
         end do

         ! angle drag-float (degrees)
         ang = editang_angat(edit_pos(2)) * 180d0 / pi
         if (edit_drag("Angle (°)##editangdrag",ang,0.2d0,2,0d0,vmax=180d0)) &
            call editang_apply(ang)
         call iw_tooltip("Angle between the three atoms (degrees)",ttshown)
         call edit_apply_button("##editang","angle")
      case (4)
         ! terminal atom combos (atoms 2 and 3 always stay fixed). The
         ! two group modes are offered independently, so imodes maps the
         ! combo slot to the move mode
         do iside = 1, 4, 3
            stropt = "Fixed"//c_null_char//"Rotate atom"//c_null_char
            nmode = 2
            imodes(1:2) = (/0,1/)
            if (w%edit_fragok(iside)) then
               stropt = stropt // "Rotate group (fix 2-3)"//c_null_char
               nmode = nmode + 1
               imodes(nmode) = 2
            end if
            ! fragok(2)/(3) = validity of the half adjacent to terminal 1/4
            if (w%edit_fragok(merge(2,3,iside == 1))) then
               stropt = stropt // "Rotate group (move 2-3)"//c_null_char
               nmode = nmode + 1
               imodes(nmode) = 3
            end if
            islot = 0
            do imode = 1, nmode
               if (imodes(imode) == w%edit_imove(iside)) islot = imode - 1
            end do
            call iw_combo_simple("Atom "//string(iside)//"##editdihmove"//string(iside),&
               stropt,islot)
            w%edit_imove(iside) = imodes(islot+1)
            call iw_tooltip("What rotates about the 2-3 axis when the dihedral changes:"//&
               " nothing (fixed), the terminal atom, the group attached to it (the"//&
               " substituents of atoms 2 and 3 stay fixed), or the whole half severed"//&
               " at the 2-3 bond (the environments of atoms 2 and 3 stay rigid)",ttshown)
         end do

         ! dihedral drag-float (degrees)
         ang = editdih_val() * 180d0 / pi
         if (edit_drag("Dihedral (°)##editdihdrag",ang,0.2d0,2,-180d0,vmax=180d0)) &
            call editdih_apply(ang)
         call iw_tooltip("Dihedral angle of the four atoms (degrees)",ttshown)
         call edit_apply_button("##editdih","dihedral")
      end select
      call panel_text("Stop editing with "//trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//&
         ", with "//trim(get_bind_keyname(BIND_CANCEL))//", with the Apply button, or by"//&
         " choosing another tool. The atoms keep the geometry they have when it stops.",&
         faded=.true.)

    end subroutine panel_edit

    ! The drag-float that sets the value of the running edit session. It
    ! is disabled while nothing is set to move, so that no value can be
    ! applied then; releasing it commits the edit. Returns whether val
    ! changed and the caller must apply it.
    function edit_drag(strid,val,speed,dec,vmin,vmax,scal) result(changed)
      character(len=*), intent(in) :: strid
      real*8, intent(inout) :: val
      real*8, intent(in) :: speed, vmin
      integer, intent(in) :: dec
      real*8, intent(in), optional :: vmax, scal
      logical :: changed

      logical :: locked, lcommit

      locked = all(w%edit_imove == 0)
      if (locked) call igBeginDisabled(.true._c_bool)
      changed = iw_dragfloat_real8(strid,x1=val,speed=speed,min=vmin,max=vmax,scale=scal,&
         decimal=dec,notlive=.true.,committed=lcommit,flags=ImGuiSliderFlags_AlwaysClamp)
      if (locked) then
         call igEndDisabled()
         changed = .false.
      elseif (lcommit) then
         call edit_commit()
      end if

    end function edit_drag

    ! The button that ends the edit session keeping the value it has
    ! reached. noun is what is being edited (distance, angle, dihedral).
    subroutine edit_apply_button(strid,noun)
      character(len=*), intent(in) :: strid, noun

      character(len=*), parameter :: natom(2:4) = (/"two  ","three","four "/)

      if (iw_button("Apply"//strid,danger=.true.)) then
         w%errmsg = ""
         call w%edit_stop()
      end if
      call iw_tooltip("Keep the current "//noun//" and release the "//&
         trim(natom(w%edit_kind))//" atoms",ttshown)

    end subroutine edit_apply_button

    ! Per-frame validity of the active edit session: the parent view
    ! must show the same system, the latched atoms must exist, dynamics
    ! must not be running, and no external geometry or bond edit may
    ! have occurred (our own edits re-stamp the heartbeat).
    function edit_session_ok() result(ok)
      logical :: ok

      ok = havesys
      if (ok) ok = (isys == w%edit_isys)
      if (ok) ok = all(w%edit_idx(1,1:w%edit_kind) >= 1) .and.&
         all(w%edit_idx(1,1:w%edit_kind) <= sys(isys)%c%ncel)
      if (ok) ok = .not.sysc(isys)%md_run
      if (ok) ok = sysc(isys)%timelastchange_geometry <= w%edit_time
      if (ok) ok = sysc(isys)%timelastchange_rebond <= w%edit_time
    end function edit_session_ok

    ! Toggle an edit session of the given kind (2 = distance, 3 =
    ! angle, 4 = dihedral): stop it if it is the active session, and
    ! start it from the measure-selected atoms otherwise.
    subroutine edit_toggle(ikind)
      integer, intent(in) :: ikind

      if (w%edit_kind == ikind) then
         call w%edit_stop()
      else
         call edit_start(ikind)
      end if
    end subroutine edit_toggle

    ! Start an edit session of the given kind (2 = distance, 3 =
    ! angle, 4 = dihedral) from the parent view's measure-selected
    ! atoms: latch the atoms and the system, compute the masked
    ! groups, and set the defaults.
    subroutine edit_start(ikind)
      integer, intent(in) :: ikind

      integer :: idx1(4), i, j
      character(len=:), allocatable :: kname

      if (.not.havesys) return
      if (.not.associated(win(iview)%sc)) return
      if (win(iview)%sc%nmsel /= ikind) return
      idx1 = 0
      idx1(1:ikind) = win(iview)%sc%msel(1,1:ikind)
      if (minval(idx1(1:ikind)) < 1 .or. maxval(idx1(1:ikind)) > sys(isys)%c%ncel) return
      do i = 1, ikind-1
         do j = i+1, ikind
            if (idx1(i) /= idx1(j)) cycle
            ! for dihedrals, atoms 2 and 3 never move, so two periodic
            ! images of the same atom are allowed there (chain crystals)
            if (ikind == 4 .and. i == 2 .and. j == 3) cycle
            w%errmsg = "Two of the selected atoms are periodic images of the same atom"
            return
         end do
      end do
      if (sysc(isys)%md_run) then
         if (ikind == 2) then
            kname = "distance"
         elseif (ikind == 3) then
            kname = "angle"
         else
            kname = "dihedral"
         end if
         w%errmsg = "Edit " // kname // " is not available while dynamics is running"
         return
      end if
      ! the session is going ahead: an armed pick mode and an edit
      ! session are mutually exclusive, so release it here, past every
      ! check that can still bail out with the mode left alone
      if (w%builder_vm /= 0) call builder_stop()
      call edit_latch(ikind)
      call edit_latch_groups(ikind)

      ! defaults: all atoms fixed except the last, which moves its
      ! group if available, else the atom; for dihedrals, prefer the
      ! group mode that also rotates the 2-3 substituents (fragok(3)
      ! is the validity of the half adjacent to terminal 4)
      w%edit_imove = 0
      if (ikind == 4 .and. w%edit_fragok(3)) then
         w%edit_imove(4) = 3
      else
         w%edit_imove(ikind) = merge(2,1,w%edit_fragok(ikind))
      end if

      w%edit_dirty = .false.
      w%edit_time = glfwGetTime()
      w%edit_kind = ikind
    end subroutine edit_start

    ! Show the latched atoms in a one-row table, one column per atom
    ! headed with the name the combos below use ("Atom 1", "Atom 2",
    ! ...). Each atom is a button colored like it is in the view. The
    ! buttons are inert: they identify the atoms, they are not controls.
    subroutine edit_atom_labels()
      integer :: is, icid
      integer(c_int) :: tflags
      real(c_float) :: rgb(3)
      logical :: havergb
      type(ImVec2) :: sz0
      character(len=:), allocatable :: lbl
      character(kind=c_char,len=:), allocatable, target :: str1

      tflags = ImGuiTableFlags_None
      tflags = ior(tflags,ImGuiTableFlags_RowBg)
      tflags = ior(tflags,ImGuiTableFlags_Borders)
      tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
      tflags = ior(tflags,ImGuiTableFlags_NoHostExtendX) ! only as wide as the columns
      str1 = "##editatomtable" // c_null_char
      sz0%x = 0._c_float
      sz0%y = 0._c_float
      if (.not.igBeginTable(c_loc(str1),int(w%edit_kind,c_int),tflags,sz0,0._c_float)) return

      do is = 1, w%edit_kind
         call iw_table_column("Atom " // string(is),id=is,flags=ImGuiTableColumnFlags_WidthFixed)
      end do
      call igTableHeadersRow()

      call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
      do is = 1, w%edit_kind
         if (.not.igTableSetColumnIndex(int(is-1,c_int))) cycle
         icid = w%edit_idx(1,is)
         ! species name, cell index and, out of the (0,0,0) cell, the lattice vector
         lbl = anchor_label(isys,int(w%edit_idx(:,is),c_int),"?",species=.true.)
         havergb = atom_view_rgb(iview,isys,atlisttype_ncel_frac,icid,rgb)
         ldum = iw_atom_button(lbl//"##editatom"//string(is),rgb,havergb,inert=.true.)
      end do
      call igEndTable()

    end subroutine edit_atom_labels

    ! Rotate the moving terminal atoms (or their groups) about the
    ! axis through center: the first latched atom rotates by
    ! -dang/nrot and the last by +dang/nrot, which opens the angle or
    ! dihedral by dang. In the dihedral "move 2-3" mode (imove = 3)
    ! the move set is the corresponding half, stored in the interior
    ! columns. The center and the axis live in the frame of the latched
    ! images, so each atom is carried to that frame (its own lattice
    ! vector within the group plus the latched vector of the seed),
    ! rotated there, and carried back: a rotation, unlike a
    ! translation, does not commute with a lattice shift.
    subroutine edit_rotate_terminals(rnew,center,axis,dang,nrot)
      real*8, intent(inout) :: rnew(:,:)
      real*8, intent(in) :: center(3), axis(3), dang
      integer, intent(in) :: nrot

      integer :: i, k, n, ia, icol
      real*8 :: rot(3,3), off(3)

      do i = 1, w%edit_kind, w%edit_kind - 1
         if (w%edit_imove(i) == 0) cycle
         rot = axisangle2mat(axis,merge(-1d0,1d0,i == 1) * dang / nrot)
         if (w%edit_imove(i) == 3) then
            icol = merge(2,3,i == 1)
            n = w%edit_nfrag(icol)
         else
            icol = i
            n = merge(w%edit_nfrag(i),1,w%edit_imove(i) == 2)
         end if
         do k = 1, n
            ia = w%edit_frag(k,icol)
            off = sys(isys)%c%x2c(dble(w%edit_fraglv(:,k,icol) + w%edit_idx(2:4,icol)))
            rnew(:,ia) = center + matmul(rot,rnew(:,ia) + off - center) - off
         end do
      end do
    end subroutine edit_rotate_terminals

    ! Latch the ikind measure-selected atoms of the parent view into
    ! edit_idx/edit_isys and clear the measurement. Stops any active
    ! session first; if its uncommitted moves rebuild the crystal, the
    ! latched lattice vectors are re-derived against the rewrapped
    ! atoms.
    subroutine edit_latch(ikind)
      integer, intent(in) :: ikind

      integer :: i, ml(4,4)
      real*8 :: rabs(3,4)

      ml = 0
      ml(:,1:ikind) = win(iview)%sc%msel(1:4,1:ikind)
      do i = 1, ikind
         rabs(:,i) = sys(isys)%c%atcel(ml(1,i))%r + sys(isys)%c%x2c(dble(ml(2:4,i)))
      end do
      call w%edit_stop()
      if (.not.sys(isys)%c%ismolecule) then
         do i = 1, ikind
            ml(2:4,i) = nint(sys(isys)%c%c2x(rabs(:,i)) - sys(isys)%c%atcel(ml(1,i))%x)
         end do
      end if
      w%edit_idx = ml
      w%edit_isys = isys
      win(iview)%sc%nmsel = 0
      win(iview)%sc%msel = 0
    end subroutine edit_latch

    ! Compute the masked groups of the ikind latched atoms: each atom
    ! severs its bonds to its neighbors in the 1-2-...-ikind chain, and
    ! its group move option is available only if the component is
    ! discrete and contains none of the other latched atoms.
    ! masked_fragment seeds each list with the atom itself, so the atom
    ! move mode is the one-atom head of the same list. For dihedrals,
    ! columns 2-3 instead hold the two halves from severing the 2-3
    ! bond: they are the move sets of the terminals' "move 2-3" modes.
    subroutine edit_latch_groups(ikind)
      integer, intent(in) :: ikind

      integer :: is, js, it, n
      logical :: okf, disc
      integer, allocatable :: ifr(:), ilv(:,:), work(:,:), worklv(:,:,:)

      allocate(work(sys(isys)%c%ncel,ikind),worklv(3,sys(isys)%c%ncel,ikind))
      w%edit_nfrag = 0
      w%edit_fragok = .false.
      do is = 1, ikind
         if (ikind == 4 .and. (is == 2 .or. is == 3)) then
            ! dihedral interior columns: the two halves obtained by
            ! severing the 2-3 bond, seeded at atoms 2 and 3. Each must
            ! contain its own terminal (1 for the 2-half, 4 for the
            ! 3-half) and no atom of the other half.
            it = merge(1,4,is == 2)
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,2),&
               w%edit_idx(1,3),n,ifr,disc,lvec=ilv)
         elseif (is > 1 .and. is < ikind) then
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is-1),n,ifr,disc,w%edit_idx(1,is),w%edit_idx(1,is+1),lvec=ilv)
         elseif (is > 1) then
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is-1),n,ifr,disc,lvec=ilv)
         else
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is+1),n,ifr,disc,lvec=ilv)
         end if
         ! valid if discrete, contains the seed's terminal (trivially
         ! true when it = is: the fragment contains its own seed), and
         ! contains no other latched atom
         okf = disc .and. any(ifr(1:n) == w%edit_idx(1,it))
         do js = 1, ikind
            if (js /= is .and. js /= it) okf = okf .and. .not.any(ifr(1:n) == w%edit_idx(1,js))
         end do
         w%edit_nfrag(is) = n
         w%edit_fragok(is) = okf
         work(1:n,is) = ifr(1:n)
         worklv(:,1:n,is) = ilv(:,1:n)
      end do
      if (allocated(w%edit_frag)) deallocate(w%edit_frag)
      if (allocated(w%edit_fraglv)) deallocate(w%edit_fraglv)
      allocate(w%edit_frag(maxval(w%edit_nfrag(1:ikind)),ikind))
      allocate(w%edit_fraglv(3,maxval(w%edit_nfrag(1:ikind)),ikind))
      w%edit_frag = 0
      w%edit_fraglv = 0
      do is = 1, ikind
         w%edit_frag(1:w%edit_nfrag(is),is) = work(1:w%edit_nfrag(is),is)
         w%edit_fraglv(:,1:w%edit_nfrag(is),is) = worklv(:,1:w%edit_nfrag(is),is)
      end do
    end subroutine edit_latch_groups

    ! Cartesian position (bohr) of the latched image of atom iside.
    function edit_pos(iside) result(r)
      integer, intent(in) :: iside
      real*8 :: r(3)

      r = sys(isys)%c%atcel(w%edit_idx(1,iside))%r +&
         sys(isys)%c%x2c(dble(w%edit_idx(2:4,iside)))
    end function edit_pos

    ! Index in the bond combo (0 = none, else bondorder(idx)) of the
    ! current bond between the two latched atoms.
    function editdist_bondidx() result(idx)
      integer :: idx
      integer :: j, k, lv(3)

      idx = 0
      lv = w%edit_idx(2:4,2) - w%edit_idx(2:4,1)
      j = sys(isys)%c%find_bond(w%edit_idx(1,1),w%edit_idx(1,2),lv)
      if (j == 0) return
      do k = 1, size(bondorder,1)
         if (sys(isys)%c%nstar(w%edit_idx(1,1))%ordcon(j) == bondorder(k)) then
            idx = k
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

      lv = w%edit_idx(2:4,2) - w%edit_idx(2:4,1)
      if (inew == 0) then
         call sysc(isys)%remove_bond(w%edit_idx(1,1),w%edit_idx(1,2),lv)
      elseif (iold == 0) then
         call sysc(isys)%add_bond(w%edit_idx(1,1),w%edit_idx(1,2),lv,bondorder(inew))
      else
         call sysc(isys)%set_bond_order(w%edit_idx(1,1),w%edit_idx(1,2),lv,bondorder(inew))
      end if
      w%edit_time = glfwGetTime()
    end subroutine editdist_setbond

    ! Move the atoms so the distance between the two latched images
    ! becomes dnew (bohr), displacing each side along the bond axis
    ! according to its move mode (in-place update; the rebuild happens
    ! at edit_commit).
    subroutine editdist_apply(dnew)
      real*8, intent(in) :: dnew

      integer :: i, k, n
      real*8 :: ra(3), rb(3), u(3), dcur, delta, dshift(2)
      real*8, allocatable :: rnew(:,:)

      ra = edit_pos(1)
      rb = edit_pos(2)
      dcur = norm2(rb - ra)
      if (dcur < eps_dzero) return
      u = (rb - ra) / dcur
      delta = dnew - dcur

      ! split the displacement among the movable sides (atom 1 moves
      ! against the bond axis, atom 2 along it)
      if (w%edit_imove(1) > 0 .and. w%edit_imove(2) > 0) then
         dshift = (/-delta/2d0, delta/2d0/)
      elseif (w%edit_imove(1) > 0) then
         dshift = (/-delta, 0d0/)
      elseif (w%edit_imove(2) > 0) then
         dshift = (/0d0, delta/)
      else
         return
      end if

      ! build the new positions and apply them in place; the moving set
      ! is the masked-fragment list (translate fragment) or its one-atom
      ! head (translate atom)
      call edit_snapshot(rnew)
      do i = 1, 2
         if (w%edit_imove(i) == 0) cycle
         n = merge(w%edit_nfrag(i),1,w%edit_imove(i) == 2)
         do k = 1, n
            rnew(:,w%edit_frag(k,i)) = rnew(:,w%edit_frag(k,i)) + dshift(i) * u
         end do
      end do
      call edit_push(rnew)
    end subroutine editdist_apply

    ! Snapshot the current atomic positions (Cartesian, bohr) into rnew.
    subroutine edit_snapshot(rnew)
      real*8, allocatable, intent(inout) :: rnew(:,:)

      integer :: i

      allocate(rnew(3,sys(isys)%c%ncel))
      do i = 1, sys(isys)%c%ncel
         rnew(:,i) = sys(isys)%c%atcel(i)%r
      end do
    end subroutine edit_snapshot

    ! Apply the moved positions in place, mark the session dirty, and
    ! notify (the rebuild happens at edit_commit).
    subroutine edit_push(rnew)
      real*8, intent(in) :: rnew(:,:)

      call sys(isys)%c%update_positions(rnew)
      w%edit_dirty = .true.
      call edit_notify()
    end subroutine edit_push

    ! Hold the camera of the parent view through the next list rebuild,
    ! post the geometry event, and re-stamp the session heartbeat
    ! (after the post, so our own edit does not end the session).
    subroutine edit_notify()

      if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
      call sysc(isys)%post_event(lastchange_geometry)
      w%edit_time = glfwGetTime()
    end subroutine edit_notify

    ! Angle (radians) at vertex position rr subtended by the two
    ! latched terminal atoms.
    function editang_angat(rr) result(a)
      real*8, intent(in) :: rr(3)
      real*8 :: a

      real*8 :: va(3), vc(3)

      va = edit_pos(1) - rr
      vc = edit_pos(3) - rr
      a = atan2(norm2(cross(va,vc)),dot_product(va,vc))
    end function editang_angat

    ! Find the translation (shift) of the central atom along the
    ! internal bisector that makes the angle equal to atgt (radians),
    ! with the terminal atoms fixed. Returns .false. if no such shift
    ! can be found.
    function editang_center_shift(atgt,shift) result(ok)
      real*8, intent(in) :: atgt
      real*8, intent(out) :: shift(3)
      logical :: ok

      integer :: i
      real*8 :: ra(3), rb(3), rc(3), va(3), vc(3), w0(3), pc(3), acur
      real*8 :: t0, t1, tm, f0, f1, fm

      ok = .false.
      shift = 0d0
      ra = edit_pos(1)
      rb = edit_pos(2)
      rc = edit_pos(3)
      va = ra - rb
      vc = rc - rb
      if (norm2(va) < eps_dzero .or. norm2(vc) < eps_dzero) return
      acur = editang_angat(rb)
      f0 = acur - atgt

      if (f0 < 0d0) then
         ! increase the angle: the internal bisector from the vertex
         ! crosses the 1-3 chord (angle-bisector theorem), where the
         ! angle tends to 180. Bracket between the vertex position and
         ! a point just short of the chord.
         pc = ra + (norm2(va) / (norm2(va) + norm2(vc))) * (rc - ra)
         t1 = norm2(pc - rb)
         if (t1 < eps_dzero) return
         w0 = (pc - rb) / t1
         t1 = t1 * (1d0 - 1d-10)
         f1 = editang_angat(rb + t1 * w0) - atgt
         if (f0 * f1 > 0d0) then
            ! the target is numerically unreachable (atgt = 180): land
            ! on the chord point, where the angle is closest to it
            shift = t1 * w0
            ok = .true.
            return
         end if
      else
         ! decrease the angle: move away from the chord along the
         ! bisector (any perpendicular to the arms if collinear); the
         ! angle goes to zero, so a doubling search brackets any
         ! reachable target
         w0 = va / norm2(va) + vc / norm2(vc)
         if (norm2(w0) < eps_dzero) then
            w0 = perpendicular(va)
         else
            w0 = w0 / norm2(w0)
         end if
         t1 = -0.05d0
         do i = 1, 60
            f1 = editang_angat(rb + t1 * w0) - atgt
            if (f0 * f1 <= 0d0) exit
            t1 = 2d0 * t1
         end do
         if (f0 * f1 > 0d0) return
      end if

      ! bisection between t0 = 0 and t1
      ok = .true.
      t0 = 0d0
      do i = 1, 80
         tm = 0.5d0 * (t0 + t1)
         fm = editang_angat(rb + tm * w0) - atgt
         if (f0 * fm > 0d0) then
            t0 = tm
            f0 = fm
         else
            t1 = tm
         end if
      end do
      shift = 0.5d0 * (t0 + t1) * w0
    end function editang_center_shift

    ! Move the atoms so the angle at the central atom becomes anew
    ! (degrees): translate the center along the bisector, or rotate
    ! the terminals about the center, according to the move modes
    ! (in-place update; the rebuild happens at edit_commit).
    subroutine editang_apply(anew)
      real*8, intent(in) :: anew

      integer :: k, n, nrot
      real*8 :: rb(3), va(3), vc(3), axis(3), acur, atgt, dang
      real*8 :: shift(3)
      real*8, allocatable :: rnew(:,:)

      rb = edit_pos(2)
      va = edit_pos(1) - rb
      vc = edit_pos(3) - rb
      if (norm2(va) < eps_dzero .or. norm2(vc) < eps_dzero) return
      acur = editang_angat(rb)
      atgt = anew * pi / 180d0
      if (abs(atgt - acur) < 1d-12) return
      if (w%edit_imove(2) > 0) then
         ! translate the central atom (or its group) along the internal
         ! bisector until the angle becomes atgt
         if (.not.editang_center_shift(atgt,shift)) return
      else
         nrot = count((/w%edit_imove(1),w%edit_imove(3)/) > 0)
         if (nrot == 0) return
         dang = atgt - acur
         ! rotation axis: perpendicular to the plane of the three atoms;
         ! any direction perpendicular to the arms if they are collinear
         axis = cross(va,vc)
         if (norm2(axis) < eps_dzero) axis = perpendicular(va)
      end if

      ! build the new positions and apply them in place; the moving set
      ! is the masked-group list or its one-atom head
      call edit_snapshot(rnew)
      if (w%edit_imove(2) > 0) then
         n = merge(w%edit_nfrag(2),1,w%edit_imove(2) == 2)
         do k = 1, n
            rnew(:,w%edit_frag(k,2)) = rnew(:,w%edit_frag(k,2)) + shift
         end do
      else
         call edit_rotate_terminals(rnew,rb,axis,dang,nrot)
      end if

      call edit_push(rnew)
    end subroutine editang_apply

    ! Signed dihedral angle (radians) of the four latched atoms about
    ! the 2-3 axis, with the same convention as the measurements
    ! (0 for degenerate, collinear input).
    function editdih_val() result(a)
      real*8 :: a

      a = sys(isys)%c%dihedral(&
         sys(isys)%c%c2x(edit_pos(1)),sys(isys)%c%c2x(edit_pos(2)),&
         sys(isys)%c%c2x(edit_pos(3)),sys(isys)%c%c2x(edit_pos(4)))
    end function editdih_val

    ! Move the atoms so the dihedral angle becomes anew (degrees):
    ! rotate the terminal atoms (or their groups) about the 2-3 axis
    ! according to the move modes (in-place update; the rebuild happens
    ! at edit_commit).
    subroutine editdih_apply(anew)
      real*8, intent(in) :: anew

      integer :: nrot
      real*8 :: r2(3), u(3), acur, dang, xn
      real*8, allocatable :: rnew(:,:)

      nrot = count((/w%edit_imove(1),w%edit_imove(4)/) > 0)
      if (nrot == 0) return

      ! the 2-3 axis; the dihedral is undefined if either arm is
      ! collinear with it
      r2 = edit_pos(2)
      u = edit_pos(3) - r2
      xn = norm2(u)
      if (xn < eps_dzero) return
      u = u / xn
      if (norm2(cross(edit_pos(1) - r2,u)) < eps_dzero) return
      if (norm2(cross(edit_pos(4) - edit_pos(3),u)) < eps_dzero) return

      acur = editdih_val()
      dang = anew * pi / 180d0 - acur
      if (abs(dang) < 1d-12) return

      ! build the new positions and apply them in place
      call edit_snapshot(rnew)
      call edit_rotate_terminals(rnew,r2,u,dang,nrot)
      call edit_push(rnew)
    end subroutine editdih_apply

    ! Commit the edit-session moves: rebuild the structure from the
    ! moved positions (keeping the current bonds) and re-latch the
    ! lattice vectors (crystal atoms may be rewrapped by the rebuild).
    subroutine edit_commit()
      integer :: i
      real*8 :: rl(3,4)

      do i = 1, w%edit_kind
         rl(:,i) = edit_pos(i)
      end do
      call sys(isys)%c%rebuild_after_move(copybonding=.true.)
      if (.not.sys(isys)%c%ismolecule) then
         do i = 1, w%edit_kind
            w%edit_idx(2:4,i) = nint(sys(isys)%c%c2x(rl(:,i)) -&
               sys(isys)%c%atcel(w%edit_idx(1,i))%x)
         end do
      end if
      w%edit_dirty = .false.
      call edit_notify()
    end subroutine edit_commit

    ! Stop the active builder mode: release the parent view (if alive)
    ! and clear the mode state.
    subroutine builder_stop()

      call bondmode_stop(w,iview,goodparent)
    end subroutine builder_stop

    ! Add-atoms mode: on empty space (icel = 0), unproject the clicked
    ! texture position xpos to the plane parallel to the screen
    ! through the scene center and add the fragment there. On an atom
    ! (icel > 0), replace it (see addatom_replace). placed is whether
    ! the structure was actually modified.
    subroutine addatom_apply(xpos,icel,placed)
      real(c_float), intent(in) :: xpos(2)
      integer, intent(in) :: icel
      logical, intent(out) :: placed

      real*8 :: x0(3), rot(3,3)
      real*8 :: tvec(3,maxaddsub), xat(3,maxaddsub+1)
      integer :: zat(maxaddsub+1), nsub, nat
      logical :: ok, used(maxaddsub)

      integer, parameter :: nobonds(2,0) = 0 ! empty attachment list

      placed = .false.
      call view_click_frame(iview,xpos,x0,rot,ok)
      if (.not.ok) return

      ! clicked on an atom: replace it instead
      if (icel > 0) then
         call addatom_replace(icel,rot,placed)
         return
      end if
      placed = .true.

      ! the fragment: central atom plus hydrogens (per-system covalent
      ! radii, same source as change_valence)
      call addatom_template(w%builder_addatom_ig,nsub,tvec)
      nat = 1
      zat(1) = w%builder_addatom_z
      xat(:,1) = x0
      used = .false.
      call addatom_cap_hydrogens(w%builder_isys,zat(1),x0,rot,nsub,tvec,used,nat,zat,xat)
      ! placed on empty space: it bonds to nothing in the system, and the
      ! system keeps the bonds it has
      call sysc(w%builder_isys)%add_atoms_fragment(nat,zat(1:nat),xat(:,1:nat),&
         newbonds=nobonds)

    end subroutine addatom_apply

    ! Add-atoms mode, click on cell atom icel: replace it with the
    ! fragment. placed is whether the structure was actually modified.
    subroutine addatom_replace(icel,rcam,placed)
      integer, intent(in) :: icel
      real*8, intent(in) :: rcam(3,3)
      logical, intent(out) :: placed

      integer :: isysl, nsub, m, ndel, nadd, k, ncon, znew, zanchor
      integer :: zat(maxaddsub+1)
      integer, allocatable :: idel(:), newbonds(:,:)
      real*8 :: tvec(3,maxaddsub), xat(3,maxaddsub+1), rot(3,3)
      real*8 :: x0(3), xtmp(3), dir(3), dnorm
      real*8, allocatable :: ufix(:,:), newbondx(:,:)
      type(substituent), allocatable :: sub(:)
      logical :: used(maxaddsub), lok

      placed = .false.
      isysl = w%builder_isys
      if (icel < 1 .or. icel > sys(isysl)%c%ncel) return
      if (.not.allocated(sys(isysl)%c%nstar)) return
      znew = w%builder_addatom_z
      call addatom_template(w%builder_addatom_ig,nsub,tvec)

      ! classify the neighbors: kept substituents vs terminal hydrogens
      ! (a mono-coordinate atom keeps its only neighbor)
      x0 = sys(isysl)%c%atcel(icel)%r
      call sys(isysl)%c%substituents(icel,ncon,sub)
      allocate(idel(ncon+1),ufix(3,max(ncon,1)),newbonds(2,max(ncon,1)),&
         newbondx(3,max(ncon,1)))
      ndel = 1
      idel(1) = icel
      m = 0
      zanchor = 0
      do k = 1, ncon
         if (ncon > 1 .and. sub(k)%z == 1 .and. sub(k)%terminal) then
            if (.not.any(idel(1:ndel) == sub(k)%id)) then
               ndel = ndel + 1
               idel(ndel) = sub(k)%id
            end if
         else
            m = m + 1
            ufix(:,m) = sys(isysl)%c%atcel(sub(k)%id)%r + sub(k)%tv
            ! the new atom takes over the bonds of the one it replaces,
            ! and makes no others: the rest of the system keeps its
            ! connectivity
            ! the image matters: in a small cell the same atom can be a
            ! neighbor through several of them
            newbonds(:,m) = (/1,sub(k)%id/)
            newbondx(:,m) = ufix(:,m)
            zanchor = sub(k)%z
         end if
      end do

      ! on a single kept substituent, re-bond the new atom at the
      ! covalent distance along the old bond direction
      if (m == 1) then
         call cov_point(isysl,ufix(:,1),zanchor,x0,znew,xtmp,lok)
         if (.not.lok) return
         x0 = xtmp
      end if

      ! enough kept substituents already: replace the element and remove the terminal hydrogens
      if (m >= nsub) then
         zat(1) = znew
         xat(:,1) = x0
         call sysc(isysl)%replace_atoms_fragment(ndel,idel(1:ndel),1,zat(1:1),xat(:,1:1),&
            newbonds=newbonds(:,1:m),newbondx=newbondx(:,1:m))
         placed = .true.
         return
      end if

      ! unit directions from the new atom to the kept substituents
      do k = 1, m
         dir = ufix(:,k) - x0
         dnorm = norm2(dir)
         if (dnorm < eps_dzero) return
         ufix(:,k) = dir / dnorm
      end do

      ! template orientation: fit the kept substituent directions, or
      ! align to the camera if there are none
      if (m > 0) then
         call addatom_fit(m,ufix,nsub,tvec,rot,used)
      else
         rot = rcam
         used = .false.
      end if

      ! the new atom plus hydrogens on the unassigned template slots
      nadd = 1
      zat(1) = znew
      xat(:,1) = x0
      call addatom_cap_hydrogens(isysl,znew,x0,rot,nsub,tvec,used,nadd,zat,xat)
      call sysc(isysl)%replace_atoms_fragment(ndel,idel(1:ndel),nadd,zat(1:nadd),xat(:,1:nadd),&
         newbonds=newbonds(:,1:m),newbondx=newbondx(:,1:m))
      placed = .true.

    end subroutine addatom_replace

    ! Append capping hydrogens at the covalent distance from the atom
    ! (element z0) at x0, on the template slots tvec not flagged in
    ! used, rotated by rot; grows the zat/xat fragment arrays from nat.
    subroutine addatom_cap_hydrogens(isysl,z0,x0,rot,nsub,tvec,used,nat,zat,xat)
      integer, intent(in) :: isysl, z0, nsub
      real*8, intent(in) :: x0(3), rot(3,3), tvec(3,maxaddsub)
      logical, intent(in) :: used(maxaddsub)
      integer, intent(inout) :: nat
      integer, intent(inout) :: zat(maxaddsub+1)
      real*8, intent(inout) :: xat(3,maxaddsub+1)

      integer :: k
      real*8 :: bl

      bl = sysc(isysl)%atmcov(z0) + sysc(isysl)%atmcov(1)
      do k = 1, nsub
         if (used(k)) cycle
         nat = nat + 1
         zat(nat) = 1
         xat(:,nat) = x0 + bl * matmul(rot,tvec(:,k))
      end do
    end subroutine addatom_cap_hydrogens

    ! Find the rotation rot and the injective assignment of the m unit
    ! directions u to the n template vertices t that best aligns them
    ! used(k) marks the template slots taken by the fitted directions.
    subroutine addatom_fit(m,u,n,t,rot,used)
      integer, intent(in) :: m
      real*8, intent(in) :: u(3,*)
      integer, intent(in) :: n
      real*8, intent(in) :: t(3,*)
      real*8, intent(out) :: rot(3,3)
      logical, intent(out) :: used(maxaddsub)

      integer :: iasg(maxaddsub), ubest(maxaddsub)
      logical :: taken(maxaddsub)
      integer :: j, ier
      real*8 :: sbest, np

      ! the identity assignment as the fallback
      rot = eye
      do j = 1, m
         iasg(j) = j
         ubest(j) = j
      end do
      call addatom_horn(m,u,t,iasg,sbest,rot,ier)
      if (ier /= 0) then
         sbest = -huge(1d0)
         rot = eye
      end if

      ! exhaustive search over the assignments if tractable
      np = 1d0
      do j = 0, m-1
         np = np * dble(n-j)
      end do
      if (np <= 1d5) then
         taken = .false.
         call addatom_fit_rec(m,u,n,t,1,iasg,taken,sbest,rot,ubest)
      end if

      used = .false.
      do j = 1, m
         used(ubest(j)) = .true.
      end do

    end subroutine addatom_fit

    ! Recursively enumerate the injective assignments of directions
    ! j..m to the free template slots, keeping the best-scoring one.
    recursive subroutine addatom_fit_rec(m,u,n,t,j,iasg,taken,sbest,rbest,ubest)
      integer, intent(in) :: m, n, j
      real*8, intent(in) :: u(3,*), t(3,*)
      integer, intent(inout) :: iasg(maxaddsub)
      logical, intent(inout) :: taken(maxaddsub)
      real*8, intent(inout) :: sbest, rbest(3,3)
      integer, intent(inout) :: ubest(maxaddsub)

      integer :: k, ier
      real*8 :: s, rot(3,3)

      do k = 1, n
         if (taken(k)) cycle
         iasg(j) = k
         if (j < m) then
            taken(k) = .true.
            call addatom_fit_rec(m,u,n,t,j+1,iasg,taken,sbest,rbest,ubest)
            taken(k) = .false.
         else
            call addatom_horn(m,u,t,iasg,s,rot,ier)
            if (ier == 0 .and. s > sbest) then
               sbest = s
               rbest = rot
               ubest(1:m) = iasg(1:m)
            end if
         end if
      end do

    end subroutine addatom_fit_rec

    ! Horn's quaternion fit (tools_math's rotation_horn) for the
    ! assignment iasg: the rotation rot maximizes
    ! sum_j u(:,j) . rot t(:,iasg(j)); the score s is m at perfect fit.
    subroutine addatom_horn(m,u,t,iasg,s,rot,ier)
      use tools_math, only: rotation_horn
      integer, intent(in) :: m
      real*8, intent(in) :: u(3,*), t(3,*)
      integer, intent(in) :: iasg(maxaddsub)
      real*8, intent(out) :: s
      real*8, intent(out) :: rot(3,3)
      integer, intent(out) :: ier

      real*8 :: tp(3,maxaddsub)
      integer :: j

      do j = 1, m
         tp(:,j) = t(:,iasg(j))
      end do
      call rotation_horn(tp(:,1:m),u(:,1:m),rot,s=s,ier=ier)

    end subroutine addatom_horn

    ! Add-atoms mode: preferred local geometry for element z, as a
    ! 0-based combo index into addgeom_names.
    function addatom_prefgeom(z) result(ig)
      integer, intent(in) :: z
      integer :: ig

      select case (z)
      case (1,3,11,19,37,55,87) ! group 1: atom
         ig = 0
      case (4,12,20,38,56,88) ! group 2: linear
         ig = 1
      case (5,13) ! B, Al: triangular
         ig = 3
      case (6,14,32) ! C, Si, Ge: tetrahedral
         ig = 6
      case (7,15,33) ! N, P, As: trigonal pyramid
         ig = 4
      case (8,16,34) ! O, S, Se: bent
         ig = 2
      case (9,17,35,53) ! F, Cl, Br, I: atom
         ig = 0
      case default ! everything else: octahedral
         ig = 10
      end select

    end function addatom_prefgeom

    ! Draw the menu entry for library fragment i and, if it is chosen, load
    ! it and arm the add-fragments mode.
    subroutine fragment_menuitem(i)
      use gui_main, only: tooltip_enabled
      integer, intent(in) :: i

      logical :: clicked

      clicked = iw_menuitem(trim(fraglib_name(i)%s))

      ! a skeletal diagram of the fragment while the entry is hovered
      if (tooltip_enabled) then
         if (igIsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByPopup)) then
            if (fragimg_ensure(i,fragslot_hover)) then
               call igBeginTooltip()
               call draw_fragment_diagram(10._c_float*igGetTextLineHeight(),fragslot_hover)
               call igEndTooltip()
            end if
         end if
      end if

      if (clicked) then
         w%errmsg = ""
         ! picking another fragment while the tool is armed only swaps
         ! the fragment: toggling would disarm the tool
         if (addfrag_load(i)) then
            w%builder_tool = vm_builder_addfragment
            if (w%builder_vm /= vm_builder_addfragment) &
               call builder_toggle(vm_builder_addfragment)
         end if
         call igCloseCurrentPopup()
      end if

    end subroutine fragment_menuitem

    ! Load fragment ifrag of the library into the builder window state.
    ! Returns .false. (and sets w%errmsg) if the fragment cannot be read
    ! or its anchor/attach atoms are not valid.
    function addfrag_load(ifrag) result(ok)
      integer, intent(in) :: ifrag
      logical :: ok

      integer :: nat, ia, it

      ! forget the fragment on deck: the panel must not keep showing the
      ! name and the diagram of the previous one next to the error
      ok = .false.
      w%builder_frag_nat = 0
      w%builder_frag_ilib = 0
      if (allocated(w%builder_frag_name)) deallocate(w%builder_frag_name)
      if (ifrag < 1 .or. ifrag > nfraglib) return

      call fraglib_geometry(ifrag,nat,w%builder_frag_z,w%builder_frag_x,ia,it,ok)
      if (.not.ok) then
         w%errmsg = "Could not read fragment " // trim(fraglib_name(ifrag)%s) //&
            " from the fragment library, or its anchor/attach atoms are invalid"
         return
      end if
      w%builder_frag_nat = nat
      w%builder_frag_ianchor = ia
      w%builder_frag_iattach = it
      w%builder_frag_radius = fraglib_radius(ifrag)
      w%builder_frag_name = trim(fraglib_name(ifrag)%s)
      w%builder_frag_ilib = ifrag

    end function addfrag_load

    ! Toggle builder mode jvm (a vm_* constant): start it on the parent
    ! view, or stop it if it is already the active mode. Starting a mode
    ! while another one is active switches to the new mode. A pick mode
    ! and an edit session cannot coexist: the first click of the new mode
    ! would kill the session anyway (the geometry event trips
    ! edit_session_ok), so end it here, where the user can see why.
    subroutine builder_toggle(jvm)
      integer, intent(in) :: jvm

      if (w%edit_kind /= 0 .and. w%builder_vm /= jvm) call w%edit_stop()
      call bondmode_toggle(w,iview,isys,jvm)
    end subroutine builder_toggle
  end subroutine draw_builder

  ! Unproject the clicked texture position xpos onto the plane
  ! parallel to the screen through the scene center: x0 is the
  ! clicked point in world coordinates and the columns of rcam are
  ! the camera axes (screen right, screen up, towards the viewer).
  subroutine view_click_frame(iview,xpos,x0,rcam,ok)
    use utils, only: invmult, mult
    use tools_math, only: cross
    integer, intent(in) :: iview
    real(c_float), intent(in) :: xpos(2)
    real*8, intent(out) :: x0(3), rcam(3,3)
    logical, intent(out) :: ok

    real(c_float) :: v0(3), vx(3), vy(3)

    ok = .false.
    if (.not.associated(win(iview)%sc)) return

    ! depth of the scene center (tworld), then unproject the click
    ! and two offset points to get the placement and camera axes
    ! (texture y points up)
    call mult(v0,win(iview)%sc%world,win(iview)%sc%scenecenter)
    call win(iview)%world_to_texpos(v0)
    vx = (/xpos(1) + 10._c_float, xpos(2), v0(3)/)
    vy = (/xpos(1), xpos(2) + 10._c_float, v0(3)/)
    v0 = (/xpos(1), xpos(2), v0(3)/)
    call win(iview)%texpos_to_world(v0)
    call win(iview)%texpos_to_world(vx)
    call win(iview)%texpos_to_world(vy)
    ! transform to world coordinates
    call invmult(v0,win(iview)%sc%world)
    call invmult(vx,win(iview)%sc%world)
    call invmult(vy,win(iview)%sc%world)
    x0 = real(v0,8)
    ! camera axes as the columns of a rotation: (x,y,z) =
    ! (screen right, screen up, toward the viewer)
    rcam(:,1) = real(vx - v0,8)
    rcam(:,1) = rcam(:,1) / norm2(rcam(:,1))
    rcam(:,2) = real(vy - v0,8)
    rcam(:,2) = rcam(:,2) - dot_product(rcam(:,2),rcam(:,1)) * rcam(:,1)
    rcam(:,2) = rcam(:,2) / norm2(rcam(:,2))
    rcam(:,3) = cross(rcam(:,1),rcam(:,2))
    ok = .true.

  end subroutine view_click_frame

  ! Add-fragments mode: place the loaded fragment. On empty space
  ! (icel = 0) the fragment goes at the clicked position, camera-aligned
  ! and capped with its placeholder atom. On an atom (icel > 0) the
  ! anchor replaces it (see addfrag_target), unless the fragment is a
  ! ligand, in which case it bonds to it (see addfrag_ligand_target).
  ! Point at the sum of the covalent radii of za (at xa) and zb, from
  ! xa towards xb (system isys). ok = the two positions are not
  ! coincident.
  subroutine cov_point(isys,xa,za,xb,zb,x,ok)
    use systems, only: sysc
    integer, intent(in) :: isys, za, zb
    real*8, intent(in) :: xa(3), xb(3)
    real*8, intent(out) :: x(3)
    logical, intent(out) :: ok

    real*8 :: dir(3), dn

    ok = .false.
    dir = xb - xa
    dn = norm2(dir)
    if (dn < eps_dzero) return
    x = xa + (sysc(isys)%atmcov(za) + sysc(isys)%atmcov(zb)) * dir / dn
    ok = .true.

  end subroutine cov_point

  subroutine frag_place(isys,iview,nat,z,x,ianchor,iattach,xdir,radius,capattach,xpos,icel,&
     nstar0,placed)
    use systems, only: sys, sysc
    use tools_math, only: cross, perpendicular
    integer, intent(in) :: isys, iview, nat, ianchor, iattach, icel
    integer, intent(in) :: z(nat)
    real*8, intent(in) :: x(3,nat), xdir(3), radius
    logical, intent(in) :: capattach
    real(c_float), intent(in) :: xpos(2)
    type(neighstar), intent(in), optional :: nstar0(nat)
    logical, intent(out), optional :: placed

    integer :: isysl, i, ia, nadd, ndel, nkept
    integer, allocatable :: idel(:), zadd(:), fmap(:), ikept(:), newbonds(:,:)
    type(neighstar), allocatable :: nstar0add(:)
    real*8, allocatable :: xadd(:,:), xkept(:,:), newbondx(:,:)
    real*8 :: x0(3), rcam(3,3), e(3,3), rot(3,3), p0(3), t1(3), dattach(3)
    logical :: ok, keepattach

    if (present(placed)) placed = .false.
    isysl = isys
    ia = ianchor
    if (nat <= 0) return
    call view_click_frame(iview,xpos,x0,rcam,ok)
    if (.not.ok) return

    ! the fragment frame and where the fragment goes
    call addfrag_frame(e)
    ! the placeholder caps the anchor on empty space, but only if it is
    ! a real atom: a dummy (the dative ligands) marks a direction only
    keepattach = .false.
    if (icel == 0 .and. iattach > 0) keepattach = (z(iattach) > 0)
    ndel = 0
    nkept = 0
    allocate(idel(1),ikept(1),xkept(3,1))
    if (icel == 0) then
       ! empty space: the attachment points to the left of the screen
       p0 = x0
       dattach = -rcam(:,1)
    elseif (radius > 0d0) then
       ! a ligand: it bonds to the clicked atom instead of replacing it
       call addfrag_ligand_target(icel,rcam,p0,dattach,ok)
       if (.not.ok) return
       nkept = 1
       ikept(1) = icel
       xkept(:,1) = sys(isysl)%c%atcel(icel)%r
    else
       call addfrag_target(icel,rcam,p0,dattach,ndel,idel,nkept,ikept,xkept,ok)
       if (.not.ok) return
    end if

    ! orient the fragment: the attachment along dattach and the plane
    ! of the fragment facing the viewer as much as it allows
    call addfrag_rotation(e,dattach,real(rcam(:,3),8),rot)
    t1 = matmul(rot,e(:,1))

    ! the fragment atoms, rotated and translated onto the anchor; the
    ! placeholder is dropped when the anchor bonds to the structure
    allocate(zadd(nat),xadd(3,nat),fmap(nat))
    nadd = 0
    do i = 1, nat
       fmap(i) = 0
       if (.not.keepattach .and. i == iattach) cycle
       nadd = nadd + 1
       fmap(i) = nadd
       zadd(nadd) = z(i)
       xadd(:,nadd) = p0 + matmul(rot,x(:,i) - x(:,ia))
    end do

    ! cap the anchor at the sum of covalent radii. Only for a library
    ! fragment, whose placeholder sits at an arbitrary distance; the atoms of
    ! a pasted fragment are real and keep the geometry they were copied with
    if (keepattach .and. capattach) &
       xadd(:,iattach) = p0 + t1 * &
       (sysc(isysl)%atmcov(z(ia)) +&
       sysc(isysl)%atmcov(z(iattach)))

    ! the internal bonding of the fragment, renumbered over the possibly
    ! dropped placeholder (its bonds die with it)
    if (present(nstar0)) &
       call nstar_subset(nstar0,fmap,nstar0add)

    ! the fragment is attached by the bonds its placement is built on and
    ! by no others: the anchor to the clicked atom (a ligand) or to the
    ! neighbors of the atom it replaces (a substituent), and nothing at
    ! all on empty space. The rest of the system keeps its bonds. With no
    ! anchor left to attach, leave newbonds unallocated (absent) and let
    ! the distance search do it
    if (fmap(ia) == 0) nkept = 0
    allocate(newbonds(2,nkept),newbondx(3,nkept))
    do i = 1, nkept
       newbonds(1,i) = fmap(ia)
       newbonds(2,i) = ikept(i)
       newbondx(:,i) = xkept(:,i)
    end do
    call sysc(isysl)%replace_atoms_fragment(ndel,idel(1:max(ndel,1)),nadd,&
       zadd(1:nadd),xadd(:,1:nadd),nstar0add,newbonds,newbondx)
    if (present(placed)) placed = .true.

  contains
    ! Add-fragments mode, click on cell atom icel: the anchor of the
    ! fragment replaces it. Returns the position of the anchor (p0), the
    ! direction from the anchor towards the atoms it bonds to (dattach),
    ! and the atoms to delete (icel plus the terminal hydrogens that make
    ! room for the fragment).
    subroutine addfrag_target(icel,rcam,p0,dattach,ndel,idel,nkept,ikept,xkept,ok)
      integer, intent(in) :: icel
      real*8, intent(in) :: rcam(3,3)
      real*8, intent(out) :: p0(3), dattach(3)
      integer, intent(out) :: ndel
      integer, allocatable, intent(inout) :: idel(:)
      integer, intent(out) :: nkept
      integer, allocatable, intent(inout) :: ikept(:)
      real*8, allocatable, intent(inout) :: xkept(:,:)
      logical, intent(out) :: ok

      integer :: isysl, ncon, k, nanch, nrem, m, zkept
      real*8 :: x0(3), dn, sumdir(3)
      type(substituent), allocatable :: sub(:)

      ok = .false.
      isysl = isys
      if (icel < 1 .or. icel > sys(isysl)%c%ncel) return
      if (.not.allocated(sys(isysl)%c%nstar)) return

      x0 = sys(isysl)%c%atcel(icel)%r
      call sys(isysl)%c%substituents(icel,ncon,sub)
      if (allocated(idel)) deallocate(idel)
      allocate(idel(ncon+1))
      ndel = 1
      idel(1) = icel
      if (allocated(ikept)) deallocate(ikept)
      if (allocated(xkept)) deallocate(xkept)
      allocate(ikept(max(ncon,1)),xkept(3,max(ncon,1)))
      nkept = 0

      ! a terminal atom is simply replaced; otherwise make room for the
      ! fragment by removing as many terminal hydrogens as connections
      ! the anchor brings with it (fewer if there are not enough)
      nanch = 0
      if (ncon > 1) nanch = addfrag_nanchor()

      ! classify the neighbors into removed hydrogens and kept neighbors
      m = 0
      nrem = 0
      zkept = 0
      sumdir = 0d0
      do k = 1, ncon
         if (sub(k)%z == 1 .and. sub(k)%terminal .and. nrem < nanch) then
            if (any(idel(1:ndel) == sub(k)%id)) cycle
            nrem = nrem + 1
            ndel = ndel + 1
            idel(ndel) = sub(k)%id
         else
            m = m + 1
            sumdir = sumdir + sub(k)%u
            zkept = sub(k)%z
            ! the anchor takes over the bonds of the atom it replaces,
            ! one per image: in a small cell the same atom can be a
            ! neighbor through more than one of them
            nkept = nkept + 1
            ikept(nkept) = sub(k)%id
            xkept(:,nkept) = sys(isysl)%c%atcel(sub(k)%id)%r + sub(k)%tv
         end if
      end do

      ! the anchor bonds towards the kept neighbors; with none left to
      ! bond to, fall back to the camera orientation
      p0 = x0
      dn = norm2(sumdir)
      if (m == 0 .or. dn < eps_dzero) then
         dattach = -rcam(:,1)
      else
         dattach = sumdir / dn
         ! on a single kept neighbor, re-bond at the sum of covalent radii
         if (m == 1) then
            call cov_point(isysl,xkept(:,1),zkept,x0,z(ianchor),p0,ok)
            if (.not.ok) return
         end if
      end if
      ok = .true.

    end subroutine addfrag_target

    ! Add-fragments mode, ligand click on cell atom icel: the whole ligand
    ! is treated as a single atom of fictitious covalent radius
    ! builder_frag_radius that bonds to icel, which is kept. Returns the
    ! position of the anchor (p0) and the direction from the anchor towards
    ! the clicked atom (dattach).
    subroutine addfrag_ligand_target(icel,rcam,p0,dattach,ok)
      integer, intent(in) :: icel
      real*8, intent(in) :: rcam(3,3)
      real*8, intent(out) :: p0(3), dattach(3)
      logical, intent(out) :: ok

      integer :: isysl, ncon, k, m, i1, i2, i3, ncand, ibest
      real*8 :: x0(3), dir(3), dn, sumdir(3), smax, sbest
      real*8 :: u(3,maxnbcand), cand(3,27)
      type(substituent), allocatable :: sub(:)

      ok = .false.
      isysl = isys
      if (icel < 1 .or. icel > sys(isysl)%c%ncel) return
      if (.not.allocated(sys(isysl)%c%nstar)) return

      ! unit directions to the neighbors of the clicked atom
      x0 = sys(isysl)%c%atcel(icel)%r
      call sys(isysl)%c%substituents(icel,ncon,sub)
      sumdir = 0d0
      m = 0
      do k = 1, ncon
         if (m >= maxnbcand) exit
         m = m + 1
         u(:,m) = sub(k)%u
         sumdir = sumdir + u(:,m)
      end do

      ! The anchor goes where there is most room, and the attachment
      ! (anchor -> clicked atom) points back the other way. Candidates are
      ! the direction opposite the sum of the bonds (exact when the shell is
      ! lopsided) plus a fixed grid in the camera frame, which still works
      ! when that sum vanishes, as it does for a linear, square planar or
      ! octahedral center. Of those, take the one whose closest existing
      ! bond is furthest away.
      if (m == 0) then
         dattach = -rcam(:,1)
      else
         ncand = 0
         dn = norm2(sumdir)
         if (dn > eps_dzero) then
            ncand = 1
            cand(:,1) = -sumdir / dn
         end if
         do i1 = -1, 1
            do i2 = -1, 1
               do i3 = -1, 1
                  if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
                  dir = real((/i1,i2,i3/),8)
                  ncand = ncand + 1
                  cand(:,ncand) = matmul(rcam,dir / norm2(dir))
               end do
            end do
         end do
         ibest = 1
         sbest = 2d0
         do k = 1, ncand
            smax = -2d0
            do i1 = 1, m
               smax = max(smax,dot_product(cand(:,k),u(:,i1)))
            end do
            if (smax < sbest) then
               sbest = smax
               ibest = k
            end if
         end do
         dattach = -cand(:,ibest)
      end if

      ! the anchor sits at the sum of the two covalent radii, and the
      ! attachment points back at the clicked atom
      p0 = x0 - (sysc(isysl)%atmcov(sys(isysl)%c%spc(sys(isysl)%c%atcel(icel)%is)%z) +&
         radius) * dattach
      ok = .true.

    end subroutine addfrag_ligand_target

    ! Number of connections the anchor brings with it when the fragment
    ! is attached: the atoms bonded to it inside the fragment, plus the
    ! attachment itself. The placeholder is counted explicitly because it
    ! may be a dummy atom, for which no covalent radius applies.
    function addfrag_nanchor() result(n)
      use global, only: bondfactor
      integer :: n

      integer :: i, ia

      ia = ianchor
      n = 1
      do i = 1, nat
         if (i == ia .or. i == iattach) cycle
         if (norm2(x(:,i) - x(:,ia)) < bondfactor *&
            (sysc(isys)%atmcov(z(i)) +&
            sysc(isys)%atmcov(z(ia)))) n = n + 1
      end do

    end function addfrag_nanchor

    ! Orthonormal frame of the loaded fragment: e1 along the attachment
    ! direction (anchor to placeholder), e3 along the direction of least
    ! spread of the fragment (the normal, for a planar fragment), and e2
    ! completing the right-handed set.
    subroutine addfrag_frame(e)
      use tools_math, only: eigsym
      real*8, intent(out) :: e(3,3)

      integer :: i, j, k, ier
      real*8 :: cov(3,3), eval(3), xd(3), dn

      ! e1: the attachment direction, from the placeholder when the fragment
      ! has one and given explicitly otherwise
      if (iattach > 0) then
         e(:,1) = x(:,iattach) - x(:,ianchor)
      else
         e(:,1) = xdir
      end if
      e(:,1) = e(:,1) / norm2(e(:,1))

      ! e3: the direction in which the fragment is least spread out
      cov = 0d0
      do i = 1, nat
         xd = x(:,i) - x(:,ianchor)
         do j = 1, 3
            do k = 1, 3
               cov(j,k) = cov(j,k) + xd(j) * xd(k)
            end do
         end do
      end do
      call eigsym(cov,3,eval,ier)
      if (ier == 0) then
         e(:,3) = cov(:,1)
      else
         e(:,3) = (/0d0,0d0,1d0/)
      end if

      ! make e3 perpendicular to e1, or any perpendicular direction if
      ! the two are parallel (a linear fragment: the twist is irrelevant)
      e(:,3) = e(:,3) - dot_product(e(:,3),e(:,1)) * e(:,1)
      dn = norm2(e(:,3))
      if (dn < eps_dzero) then
         e(:,3) = perpendicular(e(:,1))
      else
         e(:,3) = e(:,3) / dn
      end if
      e(:,2) = cross(e(:,3),e(:,1))

    end subroutine addfrag_frame

    ! Rotation mapping the fragment frame e onto the structure: the
    ! attachment direction (e1) goes along t1 and the fragment plane
    ! faces the viewer (view = the direction towards the viewer) as much
    ! as the attachment allows.
    subroutine addfrag_rotation(e,t1,view,rot)
      real*8, intent(in) :: e(3,3), t1(3), view(3)
      real*8, intent(out) :: rot(3,3)

      real*8 :: t(3,3), dn

      t(:,1) = t1 / norm2(t1)
      t(:,3) = view - dot_product(view,t(:,1)) * t(:,1)
      dn = norm2(t(:,3))
      if (dn < eps_dzero) then
         ! the attachment points at the viewer: any twist will do
         t(:,3) = perpendicular(t(:,1))
      else
         t(:,3) = t(:,3) / dn
      end if
      t(:,2) = cross(t(:,3),t(:,1))
      rot = matmul(t,transpose(e))

    end subroutine addfrag_rotation

  end subroutine frag_place

  !> Paste the clipboard fragment into system isys, as commanded by a
  !> click (or the paste key) at texture position xpos in view iview,
  !> with icel the cell atom under the cursor (0 = empty space). The
  !> fragment is placed exactly like a library fragment: on an atom the
  !> anchor replaces it, on empty space it goes at the cursor. Does
  !> nothing if the clipboard is empty.
  module subroutine paste_clipboard_fragment(isys,iview,xpos,icel,atcenter)
    use systems, only: sysclip, ok_system, sys_init
    use utils, only: mult
    use tools_math, only: m_x2c_from_cellpar
    integer, intent(in) :: isys, iview, icel
    real(c_float), intent(in) :: xpos(2)
    logical, intent(in), optional :: atcenter

    integer :: i, nat
    integer, allocatable :: z(:), imap(:)
    real*8, allocatable :: x(:,:)
    real*8 :: x2c(3,3)
    real(c_float) :: xpos_(2), v0(3)
    type(neighstar), allocatable :: nstar0(:)

    ! consistency checks
    if (.not.sysclip%isfilled) return
    if (.not.ok_system(isys,sys_init)) return
    nat = sysclip%seed%nat
    if (nat <= 0) return
    if (sysclip%ianchor < 1 .or. sysclip%ianchor > nat) return

    ! the clipboard atoms in Cartesian (bohr): a fragment copied from a
    ! crystal is stored in fractional coordinates
    allocate(z(nat),x(3,nat))
    x2c = 0d0
    if (sysclip%seed%useabr == 1) then
       x2c = m_x2c_from_cellpar(sysclip%seed%aa,sysclip%seed%bb)
    elseif (sysclip%seed%useabr == 2) then
       x2c = sysclip%seed%m_x2c
    end if
    do i = 1, nat
       z(i) = sysclip%seed%spc(sysclip%seed%is(i))%z
       if (sysclip%seed%useabr > 0) then
          x(:,i) = matmul(x2c,sysclip%seed%x(:,i))
       else
          x(:,i) = sysclip%seed%x(:,i)
       end if
    end do

    ! when commanded from the menu there is no cursor to place it at, so use
    ! the middle of the view (the texture position of the scene center)
    xpos_ = xpos
    if (present(atcenter)) then
       if (atcenter) then
          if (.not.associated(win(iview)%sc)) return
          call mult(v0,win(iview)%sc%world,win(iview)%sc%scenecenter)
          call win(iview)%world_to_texpos(v0)
          xpos_ = (/v0(1),v0(2)/)
       end if
    end if

    ! The internal bonding of the copied fragment, so the paste restores
    ! it instead of recomputing it from the distances. The pasted atoms
    ! form a single piece: bonds to other periodic images are dropped.
    if (sysclip%seed%havebonds) then
       allocate(imap(nat))
       do i = 1, nat
          imap(i) = i
       end do
       call nstar_subset(sysclip%seed%nstar,imap,nstar0,droppbc=.true.)
    end if

    ! a pasted fragment is always a substituent, never a ligand
    call frag_place(isys,iview,nat,z,x,sysclip%ianchor,sysclip%iattach,&
       sysclip%xdir,0d0,.false.,xpos_,icel,nstar0)

  end subroutine paste_clipboard_fragment

  !> Stop the active edit session (distance, angle, or dihedral) of a builder
  !> window: rebuild the structure first if in-place moves were applied
  !> but not committed, and release the latched atoms (the start
  !> routines re-seed the rest of the session state).
  module subroutine edit_stop(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_geometry
    class(window), intent(inout) :: w

    integer :: iview

    if (w%edit_dirty) then
       if (ok_system(w%edit_isys,sys_init)) then
          call sys(w%edit_isys)%c%rebuild_after_move(copybonding=.true.)
          iview = w%idparent
          if (iview >= 1 .and. iview <= nwin) then
             if (win(iview)%isinit) then
                if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
             end if
          end if
          call sysc(w%edit_isys)%post_event(lastchange_geometry)
       end if
       w%edit_dirty = .false.
    end if
    w%edit_kind = 0
    w%edit_isys = 0
    w%edit_idx = 0

  end subroutine edit_stop

  ! Add-atoms mode: number of substituents n and unit vectors v(:,1:n)
  ! for local geometry index ig (0-based, see addgeom_names). Local
  ! frame: x = screen right, y = screen up, z = toward the viewer.
  subroutine addatom_template(ig,n,v)
    integer, intent(in) :: ig
    integer, intent(out) :: n
    real*8, intent(out) :: v(3,maxaddsub)

    integer :: k

    real*8, parameter :: ath = 109.4712d0 ! tetrahedral angle (degrees)

    n = 0
    v = 0d0
    select case (ig)
    case (1) ! linear
       n = 2
       call addatom_sphvec(90d0,0d0,v(:,1))
       call addatom_sphvec(90d0,180d0,v(:,2))
    case (2) ! bent
       n = 2
       call addatom_sphvec(ath/2d0,0d0,v(:,1))
       call addatom_sphvec(ath/2d0,180d0,v(:,2))
    case (3) ! triangular
       n = 3
       call addatom_sphvec(0d0,0d0,v(:,1))
       call addatom_sphvec(120d0,0d0,v(:,2))
       call addatom_sphvec(120d0,180d0,v(:,3))
    case (4) ! trigonal pyramid
       n = 3
       do k = 1, 3
          call addatom_sphvec(ath,90d0+120d0*(k-1),v(:,k))
       end do
    case (5) ! T-shape
       n = 3
       call addatom_sphvec(0d0,0d0,v(:,1))
       call addatom_sphvec(180d0,0d0,v(:,2))
       call addatom_sphvec(90d0,180d0,v(:,3))
    case (6) ! tetrahedral
       n = 4
       call addatom_sphvec(0d0,0d0,v(:,1))
       do k = 1, 3
          call addatom_sphvec(ath,90d0+120d0*(k-1),v(:,1+k))
       end do
    case (7) ! square planar
       n = 4
       do k = 1, 4
          call addatom_sphvec(90d0,45d0+90d0*(k-1),v(:,k))
       end do
    case (8) ! trigonal bipyramid
       n = 5
       call addatom_sphvec(0d0,0d0,v(:,1))
       call addatom_sphvec(180d0,0d0,v(:,2))
       do k = 1, 3
          call addatom_sphvec(90d0,90d0+120d0*(k-1),v(:,2+k))
       end do
    case (9) ! square pyramid
       n = 5
       call addatom_sphvec(0d0,0d0,v(:,1))
       do k = 1, 4
          call addatom_sphvec(100d0,45d0+90d0*(k-1),v(:,1+k))
       end do
    case (10) ! octahedral
       n = 6
       call addatom_sphvec(0d0,0d0,v(:,1))
       call addatom_sphvec(180d0,0d0,v(:,2))
       do k = 1, 4
          call addatom_sphvec(90d0,90d0*(k-1),v(:,2+k))
       end do
    case (11) ! trigonal antiprism
       n = 6
       do k = 1, 3
          call addatom_sphvec(60d0,30d0+120d0*(k-1),v(:,k))
          call addatom_sphvec(120d0,90d0+120d0*(k-1),v(:,3+k))
       end do
    case (12) ! pentagonal bipyramid
       n = 7
       call addatom_sphvec(0d0,0d0,v(:,1))
       call addatom_sphvec(180d0,0d0,v(:,2))
       do k = 1, 5
          call addatom_sphvec(90d0,90d0+72d0*(k-1),v(:,2+k))
       end do
    case (13) ! capped octahedron
       n = 7
       call addatom_sphvec(0d0,0d0,v(:,1))
       do k = 1, 3
          call addatom_sphvec(75d0,30d0+120d0*(k-1),v(:,1+k))
          call addatom_sphvec(130d0,90d0+120d0*(k-1),v(:,4+k))
       end do
    case (14) ! square antiprism
       n = 8
       do k = 1, 4
          call addatom_sphvec(59d0,90d0*(k-1),v(:,k))
          call addatom_sphvec(121d0,45d0+90d0*(k-1),v(:,4+k))
       end do
    case (15) ! trigonal dodecahedron
       n = 8
       call addatom_sphvec(36.9d0,0d0,v(:,1))
       call addatom_sphvec(36.9d0,180d0,v(:,2))
       call addatom_sphvec(143.1d0,90d0,v(:,3))
       call addatom_sphvec(143.1d0,270d0,v(:,4))
       call addatom_sphvec(69.5d0,90d0,v(:,5))
       call addatom_sphvec(69.5d0,270d0,v(:,6))
       call addatom_sphvec(110.5d0,0d0,v(:,7))
       call addatom_sphvec(110.5d0,180d0,v(:,8))
    case (16) ! tricapped trigonal prism
       n = 9
       do k = 1, 3
          call addatom_sphvec(48d0,30d0+120d0*(k-1),v(:,k))
          call addatom_sphvec(132d0,30d0+120d0*(k-1),v(:,3+k))
          call addatom_sphvec(90d0,90d0+120d0*(k-1),v(:,6+k))
       end do
    case (17) ! capped square antiprism
       n = 9
       call addatom_sphvec(0d0,0d0,v(:,1))
       do k = 1, 4
          call addatom_sphvec(68d0,90d0*(k-1),v(:,1+k))
          call addatom_sphvec(118d0,45d0+90d0*(k-1),v(:,5+k))
       end do
    case (18) ! pentagonal prism
       n = 10
       do k = 1, 5
          call addatom_sphvec(55d0,90d0+72d0*(k-1),v(:,k))
          call addatom_sphvec(125d0,90d0+72d0*(k-1),v(:,5+k))
       end do
    end select

  end subroutine addatom_template

  ! Unit vector at polar angle th (degrees, from screen up) and
  ! azimuth ph (degrees, from screen right towards the viewer).
  subroutine addatom_sphvec(th,ph,v)
    use param, only: pi
    real*8, intent(in) :: th, ph
    real*8, intent(out) :: v(3)

    real*8 :: t, p

    t = th * pi / 180d0
    p = ph * pi / 180d0
    v = (/sin(t)*cos(p), cos(t), sin(t)*sin(p)/)

  end subroutine addatom_sphvec

  !> Draw the add-atoms local-geometry selector: a grid of square
  !> pictograms, one per geometry. ig is the current selection (0-based
  !> index), updated when an icon is clicked. sidefac is the side of a
  !> pictogram, in text line heights. Hovering an icon shows the
  !> geometry name.
  subroutine draw_addatom_geom_grid(ig,sidefac)
    use gui_main, only: g
    use utils, only: iw_tooltip
    integer, intent(inout) :: ig
    real(c_float), intent(in) :: sidefac

    integer :: i, ncol
    logical :: hovered
    type(ImVec2) :: sz, p0, p1, szavail
    type(c_ptr) :: dl
    integer(c_int) :: col
    real(c_float) :: side
    character(len=:,kind=c_char), allocatable, target :: str1

    side = sidefac * igGetTextLineHeightWithSpacing()

    ! as many icons per row as fit in the available width
    call igGetContentRegionAvail(szavail)
    ncol = int((szavail%x + g%Style%ItemInnerSpacing%x) / (side + g%Style%ItemInnerSpacing%x))
    ncol = max(min(ncol,naddgeom),4)
    sz%x = side
    sz%y = side
    str1 = "##geomicon" // c_null_char
    dl = igGetWindowDrawList()
    do i = 0, naddgeom-1
       if (mod(i,ncol) /= 0) call igSameLine(0._c_float,g%Style%ItemInnerSpacing%x)
       call igPushID_Int(int(i,c_int))
       if (igInvisibleButton(c_loc(str1),sz,ImGuiButtonFlags_None)) ig = i
       hovered = igIsItemHovered(ImGuiHoveredFlags_None)
       call igGetItemRectMin(p0)
       p1%x = p0%x + side
       p1%y = p0%y + side

       ! background, frame, and the pictogram
       if (ig == i) then
          col = igGetColorU32_Col(ImGuiCol_ButtonActive,0.6_c_float)
          call ImDrawList_AddRectFilled(dl,p0,p1,col,0._c_float,0_c_int)
       elseif (hovered) then
          col = igGetColorU32_Col(ImGuiCol_ButtonHovered,0.6_c_float)
          call ImDrawList_AddRectFilled(dl,p0,p1,col,0._c_float,0_c_int)
       end if
       col = igGetColorU32_Col(merge(ImGuiCol_Text,ImGuiCol_Border,ig == i),1._c_float)
       call ImDrawList_AddRect(dl,p0,p1,col,0._c_float,0_c_int,1._c_float)
       call addatom_geom_paint(i,p0,side,dl)

       ! the geometry name in a tooltip
       if (hovered) call iw_tooltip(trim(addgeom_names(i+1)))
       call igPopID()
    end do

  end subroutine draw_addatom_geom_grid

  !> Paint the pictogram for add-atoms local geometry ig (0-based)
  !> into the square with top-left screen position p0 and side length
  !> side, using draw list dl. If iz is present and positive, the
  !> central atom is the symbol of element iz in the element color
  !> instead of a filled circle. If bgrgb is present, choose the
  !> colors to contrast with a backdrop of that color (otherwise the
  !> ImGui text color is used).
  module subroutine addatom_geom_paint(ig,p0,side,dl,iz,bgrgb)
    use tools_io, only: nameguess
    use param, only: pi, jmlcol
    integer, intent(in) :: ig
    type(ImVec2), intent(in) :: p0
    real(c_float), intent(in) :: side
    type(c_ptr), intent(in) :: dl
    integer, intent(in), optional :: iz
    real(c_float), intent(in), optional :: bgrgb(3)

    logical :: dosym, dark
    integer :: nb, k, j
    integer :: sty(maxaddsub)
    real*8 :: ang(maxaddsub), lfac(maxaddsub), a
    type(ImVec2) :: cen, pa, pb, q1, q2, tsz
    type(ImVec4) :: cv
    integer(c_int) :: col
    real(c_float) :: ux, uy, blen, f, hw, rstart, rmax, sluma
    character(len=:,kind=c_char), allocatable, target :: str1

    real(c_float), parameter :: r0rat = 0.14_c_float ! central circle radius (fraction of side)
    real(c_float), parameter :: blrat = 0.40_c_float ! bond length (fraction of side)
    real(c_float), parameter :: wedgerat = 0.09_c_float ! wedge/hash half-width at the far end (fraction of side)
    integer, parameter :: nhash = 4 ! number of strokes in a hashed wedge
    real(c_float), parameter :: lumamin = 0.45_c_float ! minimum symbol luminance on a dark backdrop
    real(c_float), parameter :: lumamax = 0.60_c_float ! maximum symbol luminance on a light backdrop

    call addatom_diagram(ig,nb,ang,sty,lfac)
    cen%x = p0%x + 0.5_c_float * side
    cen%y = p0%y + 0.5_c_float * side

    ! bond color: black or white against the given backdrop, or the ImGui text color
    if (present(bgrgb)) then
       cv%x = 0._c_float
       cv%y = 0._c_float
       cv%z = 0._c_float
       cv%w = 1._c_float
       if (luma(bgrgb(1),bgrgb(2),bgrgb(3)) <= 0.5_c_float) then
          cv%x = 1._c_float
          cv%y = 1._c_float
          cv%z = 1._c_float
       end if
       col = igGetColorU32_Vec4(cv)
    else
       col = igGetColorU32_Col(ImGuiCol_Text,1._c_float)
    end if

    ! the bonds start at the edge of the central circle or, if the
    ! element symbol is shown, outside its text box
    dosym = .false.
    if (present(iz)) dosym = (iz > 0)
    rstart = r0rat * side
    if (dosym) then
       str1 = trim(nameguess(iz,.true.)) // c_null_char
       call igCalcTextSize(tsz,c_loc(str1),c_null_ptr,.false._c_bool,-1._c_float)
       rstart = max(rstart,0.55_c_float * sqrt(tsz%x*tsz%x + tsz%y*tsz%y))
    end if

    ! the bonds must end inside the box, wedges included. The symbol is
    ! drawn at the font size, so in a small box it can push the start of
    ! the bonds so far out that the tips would stick out of it
    rmax = (0.5_c_float - wedgerat) * side
    rstart = min(rstart,rmax - 0.12_c_float * side)

    do k = 1, nb
       ! bond direction on screen (angles counterclockwise from +x, screen y points down)
       a = ang(k) * pi / 180d0
       ux = real(cos(a),c_float)
       uy = -real(sin(a),c_float)
       blen = min(rstart + (blrat * real(lfac(k),c_float) - r0rat) * side,rmax)
       pa%x = cen%x + rstart * ux
       pa%y = cen%y + rstart * uy
       pb%x = cen%x + blen * ux
       pb%y = cen%y + blen * uy
       if (sty(k) == sty_wedge) then
          ! towards the viewer: filled wedge, wide at the far end
          q1%x = pb%x - wedgerat * side * uy
          q1%y = pb%y + wedgerat * side * ux
          q2%x = pb%x + wedgerat * side * uy
          q2%y = pb%y - wedgerat * side * ux
          call ImDrawList_AddTriangleFilled(dl,pa,q1,q2,col)
       elseif (sty(k) == sty_hash) then
          ! away from the viewer: hashed wedge, strokes perpendicular
          ! to the bond that widen with the distance
          do j = 1, nhash
             f = 0.30_c_float + 0.70_c_float * real(j-1,c_float) / real(nhash-1,c_float)
             hw = f * wedgerat * side
             q1%x = pa%x + f * (pb%x - pa%x) - hw * uy
             q1%y = pa%y + f * (pb%y - pa%y) + hw * ux
             q2%x = pa%x + f * (pb%x - pa%x) + hw * uy
             q2%y = pa%y + f * (pb%y - pa%y) - hw * ux
             call ImDrawList_AddLine(dl,q1,q2,col,max(1._c_float,0.045_c_float*side))
          end do
       else
          ! in the plane of the diagram: plain line
          call ImDrawList_AddLine(dl,pa,pb,col,max(1.5_c_float,0.06_c_float*side))
       end if
    end do

    ! the central atom: element symbol in the element color, or a filled circle
    if (dosym) then
       cv%x = real(jmlcol(1,iz),c_float) / 255._c_float
       cv%y = real(jmlcol(2,iz),c_float) / 255._c_float
       cv%z = real(jmlcol(3,iz),c_float) / 255._c_float
       cv%w = 1._c_float

       ! clamp the luminance so the symbol contrasts with the
       ! backdrop: darken light colors on a light backdrop, lighten
       ! dark colors on a dark one
       sluma = luma(cv%x,cv%y,cv%z)
       if (present(bgrgb)) then
          dark = (luma(bgrgb(1),bgrgb(2),bgrgb(3)) <= 0.5_c_float)
       else
          dark = .true.
       end if
       if (dark) then
          if (sluma < lumamin) then
             cv%x = min(cv%x + lumamin - sluma,1._c_float)
             cv%y = min(cv%y + lumamin - sluma,1._c_float)
             cv%z = min(cv%z + lumamin - sluma,1._c_float)
          end if
       else
          if (sluma > lumamax) then
             cv%x = cv%x * lumamax / sluma
             cv%y = cv%y * lumamax / sluma
             cv%z = cv%z * lumamax / sluma
          end if
       end if

       q1%x = cen%x - 0.5_c_float * tsz%x
       q1%y = cen%y - 0.5_c_float * tsz%y
       call ImDrawList_AddText_Vec2(dl,q1,igGetColorU32_Vec4(cv),c_loc(str1),c_null_ptr)
    else
       call ImDrawList_AddCircleFilled(dl,cen,r0rat*side,col,0_c_int)
    end if

  contains
    ! relative luminance of an rgb color (Rec. 709 weights)
    function luma(r,g,b)
      real(c_float), intent(in) :: r, g, b
      real(c_float) :: luma
      luma = 0.2126_c_float * r + 0.7152_c_float * g + 0.0722_c_float * b
    end function luma
  end subroutine addatom_geom_paint

  ! Return the 2D diagram for local atom geometry ig (0-based): nb
  ! bonds with screen angles ang (degrees, counterclockwise from +x),
  ! styles sty (plain line in the plane of the diagram, filled wedge
  ! towards the viewer, or hashed wedge away from the viewer), and
  ! bond length factors lfac.
  subroutine addatom_diagram(ig,nb,ang,sty,lfac)
    integer, intent(in) :: ig
    integer, intent(out) :: nb
    real*8, intent(out) :: ang(maxaddsub), lfac(maxaddsub)
    integer, intent(out) :: sty(maxaddsub)

    nb = 0
    ang = 0d0
    sty = sty_plain
    lfac = 1d0
    select case(ig)
    case(1) ! linear
       nb = 2
       ang(1:2) = (/0d0,180d0/)
    case(2) ! bent
       nb = 2
       ang(1:2) = (/55d0,125d0/)
    case(3) ! triangular
       nb = 3
       ang(1:3) = (/90d0,210d0,330d0/)
    case(4) ! trigonal pyramid
       nb = 3
       ang(1:3) = (/210d0,270d0,330d0/)
       sty(1:3) = (/sty_plain,sty_hash,sty_wedge/)
    case(5) ! T-shape
       nb = 3
       ang(1:3) = (/0d0,180d0,270d0/)
    case(6) ! tetrahedral
       nb = 4
       ang(1:4) = (/90d0,210d0,270d0,330d0/)
       sty(1:4) = (/sty_plain,sty_plain,sty_hash,sty_wedge/)
    case(7) ! square planar
       nb = 4
       ang(1:4) = (/0d0,180d0,135d0,315d0/)
       sty(1:4) = (/sty_plain,sty_plain,sty_hash,sty_wedge/)
    case(8) ! trigonal bipyramid
       nb = 5
       ang(1:5) = (/90d0,270d0,180d0,20d0,340d0/)
       sty(1:5) = (/sty_plain,sty_plain,sty_plain,sty_hash,sty_wedge/)
    case(9) ! square pyramid
       nb = 5
       ang(1:5) = (/90d0,30d0,150d0,210d0,330d0/)
       sty(1:5) = (/sty_plain,sty_hash,sty_hash,sty_wedge,sty_wedge/)
    case(10) ! octahedral
       nb = 6
       ang(1:6) = (/90d0,270d0,30d0,150d0,210d0,330d0/)
       sty(1:6) = (/sty_plain,sty_plain,sty_hash,sty_hash,sty_wedge,sty_wedge/)
    case(11) ! trigonal antiprism
       nb = 6
       ang(1:6) = (/50d0,90d0,130d0,230d0,270d0,310d0/)
       sty(1:3) = sty_hash
       sty(4:6) = sty_wedge
    case(12) ! pentagonal bipyramid
       nb = 7
       ang(1:7) = (/90d0,270d0,270d0,0d0,180d0,216d0,324d0/)
       sty(1:7) = (/sty_plain,sty_plain,sty_hash,sty_plain,sty_plain,sty_wedge,sty_wedge/)
       lfac(2:3) = (/1.15d0,0.7d0/)
    case(13) ! capped octahedron
       nb = 7
       ang(1:7) = (/90d0,90d0,25d0,155d0,210d0,270d0,330d0/)
       sty(1:7) = (/sty_plain,sty_hash,sty_hash,sty_hash,sty_wedge,sty_wedge,sty_wedge/)
       lfac(1:2) = (/1.15d0,0.7d0/)
    case(14) ! square antiprism
       nb = 8
       ang(1:8) = (/22d0,68d0,112d0,158d0,202d0,248d0,292d0,338d0/)
       sty(1:4) = sty_hash
       sty(5:8) = sty_wedge
    case(15) ! trigonal dodecahedron
       nb = 8
       ang(1:8) = (/15d0,75d0,105d0,165d0,195d0,255d0,285d0,345d0/)
       sty(1:4) = sty_hash
       sty(5:8) = sty_wedge
    case(16) ! tricapped trigonal prism
       nb = 9
       ang(1:9) = (/0d0,180d0,90d0,90d0,40d0,140d0,220d0,270d0,320d0/)
       sty(1:9) = (/sty_plain,sty_plain,sty_plain,sty_hash,sty_hash,sty_hash,&
          sty_wedge,sty_wedge,sty_wedge/)
       lfac(3:4) = (/1.15d0,0.7d0/)
    case(17) ! capped square antiprism
       nb = 9
       ang(1:9) = (/90d0,20d0,65d0,115d0,160d0,205d0,250d0,295d0,340d0/)
       sty(2:5) = sty_hash
       sty(6:9) = sty_wedge
       lfac(1) = 1.15d0
    case(18) ! pentagonal prism
       nb = 10
       ang(1:10) = (/18d0,54d0,90d0,126d0,162d0,198d0,234d0,270d0,306d0,342d0/)
       sty(1:5) = sty_hash
       sty(6:10) = sty_wedge
    end select

  end subroutine addatom_diagram

  ! Read the coordinates of library fragment ifrag. z and x (Cartesian,
  ! bohr) are allocated here; ia and it are the anchor and the attachment
  ! placeholder, validated against the atom count.
  subroutine fraglib_geometry(ifrag,nat,z,x,ia,it,ok)
    use crystalseedmod, only: crystalseed
    use global, only: flib_file
    integer, intent(in) :: ifrag
    integer, intent(out) :: nat, ia, it
    integer, allocatable, intent(inout) :: z(:)
    real*8, allocatable, intent(inout) :: x(:,:)
    logical, intent(out) :: ok

    type(crystalseed) :: seed
    integer :: i

    ok = .false.
    nat = 0
    if (ifrag < 1 .or. ifrag > nfraglib) return

    call seed%read_library(trim(fraglib_name(ifrag)%s),ok,file=flib_file)
    if (.not.ok .or. seed%nat <= 0 .or. seed%nspc <= 0) then
       ok = .false.
       return
    end if

    ia = fraglib_anchor(ifrag)
    it = fraglib_attach(ifrag)
    ok = (ia >= 1 .and. ia <= seed%nat .and. it >= 1 .and. it <= seed%nat .and. ia /= it)
    if (.not.ok) return

    if (allocated(z)) deallocate(z)
    if (allocated(x)) deallocate(x)
    allocate(z(seed%nat),x(3,seed%nat))
    do i = 1, seed%nat
       z(i) = seed%spc(seed%is(i))%z
       x(:,i) = seed%x(:,i)
    end do
    nat = seed%nat

  end subroutine fraglib_geometry

  ! Draw the structural diagram held in texture slot islot as an inline
  ! image of the given side length. The diagrams are pre-rendered
  ! (dat/assets/fragments) and loaded on demand by fragimg_ensure.
  subroutine draw_fragment_diagram(side,islot)
    use interfaces_opengl3
    real(c_float), intent(in) :: side
    integer, intent(in) :: islot

    type(ImVec2) :: sz, uv0, uv1
    type(ImVec4) :: tint, nobord

    if (fragimg_tex(islot) == 0) return
    sz%x = side
    sz%y = side
    uv0%x = 0._c_float
    uv0%y = 0._c_float
    uv1%x = 1._c_float
    uv1%y = 1._c_float
    tint = ImVec4(1._c_float,1._c_float,1._c_float,1._c_float)
    nobord = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)
    call igImage(int(fragimg_tex(islot),c_intptr_t),sz,uv0,uv1,tint,nobord)

  end subroutine draw_fragment_diagram

  ! Load the diagram of library fragment ifrag into texture slot islot.
  ! Returns .false. if the fragment has no image, so the caller can skip
  ! the tooltip instead of opening an empty one.
  function fragimg_ensure(ifrag,islot)
    use interfaces_stb, only: stbi_load, stbi_image_free
    use interfaces_opengl3
    use global, only: critic_home
    use param, only: dirsep
    integer, intent(in) :: ifrag, islot
    logical :: fragimg_ensure

    integer(c_int) :: iwidth, iheight, ichannel
    type(c_ptr) :: pixels
    character(kind=c_char,len=:), allocatable, target :: file

    fragimg_ensure = .false.
    if (ifrag < 1 .or. ifrag > nfraglib) return
    if (fragimg_idx(islot) == ifrag) then
       fragimg_ensure = (fragimg_tex(islot) /= 0)
       return
    end if

    if (fragimg_tex(islot) /= 0) call glDeleteTextures(1,c_loc(fragimg_tex(islot)))
    fragimg_tex(islot) = 0
    fragimg_idx(islot) = ifrag
    file = trim(critic_home) // dirsep // "assets" // dirsep // "fragments" //&
       dirsep // trim(fraglib_name(ifrag)%s) // ".png" // c_null_char
    pixels = stbi_load(c_loc(file),iwidth,iheight,ichannel,4)
    if (.not.c_associated(pixels)) return

    call glGenTextures(1,c_loc(fragimg_tex(islot)))
    call glBindTexture(GL_TEXTURE_2D,fragimg_tex(islot))
    call glPixelStorei(GL_UNPACK_ALIGNMENT,1)
    call glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,iwidth,iheight,0,GL_RGBA,&
       GL_UNSIGNED_BYTE,pixels)
    call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR)
    call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR)
    call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE)
    call glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE)
    call glBindTexture(GL_TEXTURE_2D,0)
    call stbi_image_free(pixels)
    fragimg_ensure = .true.

  end function fragimg_ensure

  ! Whether fragment j is the first in the library with its category and,
  ! if lsub, also with its subcategory: marks where a submenu is opened.
  function fraglib_isfirst(j,lsub)
    integer, intent(in) :: j
    logical, intent(in) :: lsub
    logical :: fraglib_isfirst

    integer :: i

    fraglib_isfirst = .false.
    do i = 1, j-1
       if (fraglib_cat(i)%s /= fraglib_cat(j)%s) cycle
       if (lsub) then
          if (fraglib_sub(i)%s /= fraglib_sub(j)%s) cycle
       end if
       return
    end do
    fraglib_isfirst = .true.

  end function fraglib_isfirst

  ! Read the list of fragments in the fragment library (flib_file) the
  ! first time it is needed: the name of each fragment plus the indices
  ! of its anchor and attachment-placeholder atoms. The format is the
  ! same as in the crystal and molecule libraries, with the anchor and
  ! attach keywords given after the endmolecule line (see the file).
  subroutine fraglib_ensure()
    use tools_io, only: fopen_read, fclose, getline, lgetword, getword, equal,&
       isinteger, isreal
    use types, only: realloc
    use global, only: flib_file
    use param, only: bohrtoa

    integer :: lu, lp, idum
    real*8 :: rdum
    character(len=:), allocatable :: line, word, name
    logical :: ok

    if (fraglib_init) return
    fraglib_init = .true.
    nfraglib = 0
    if (.not.allocated(flib_file)) return
    inquire(file=flib_file,exist=ok)
    if (.not.ok) return
    lu = fopen_read(flib_file,abspath0=.true.)
    if (lu < 0) return

    allocate(fraglib_name(10),fraglib_anchor(10),fraglib_attach(10),fraglib_cat(10),&
       fraglib_sub(10),fraglib_radius(10))
    do while (getline(lu,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"structure")) then
          ! a new fragment: default to the first atom as the anchor and
          ! no placeholder until the attach keyword is found
          name = getword(line,lp)
          if (len_trim(name) == 0) cycle
          nfraglib = nfraglib + 1
          if (nfraglib > size(fraglib_name,1)) then
             call realloc(fraglib_name,2*nfraglib)
             call realloc(fraglib_anchor,2*nfraglib)
             call realloc(fraglib_attach,2*nfraglib)
             call realloc(fraglib_cat,2*nfraglib)
             call realloc(fraglib_sub,2*nfraglib)
             call realloc(fraglib_radius,2*nfraglib)
          end if
          fraglib_name(nfraglib)%s = name
          fraglib_anchor(nfraglib) = 1
          fraglib_attach(nfraglib) = 0
          fraglib_cat(nfraglib)%s = ""
          fraglib_sub(nfraglib)%s = ""
          fraglib_radius(nfraglib) = 0d0
       elseif (nfraglib > 0) then
          if (equal(word,"anchor")) then
             if (isinteger(idum,line,lp)) fraglib_anchor(nfraglib) = idum
          elseif (equal(word,"attach")) then
             if (isinteger(idum,line,lp)) fraglib_attach(nfraglib) = idum
          elseif (equal(word,"radius")) then
             ! a fictitious covalent radius marks a ligand: the whole fragment
             ! is placed like a single atom of that radius, bonded to the
             ! clicked metal instead of replacing it (the file gives it in
             ! angstrom, everything downstream is in bohr)
             if (isreal(rdum,line,lp)) fraglib_radius(nfraglib) = rdum / bohrtoa
          elseif (equal(word,"category")) then
             ! the category is free text (it labels a submenu); a category of
             ! the form "parent/child" nests the fragment one level deeper
             ! split only when both halves are non-empty: an empty menu
             ! label is fatal in imgui, and this file is user-editable
             name = trim(adjustl(line(lp:)))
             idum = index(name,"/")
             if (idum > 1 .and. len_trim(name(idum+1:)) > 0) then
                fraglib_cat(nfraglib)%s = trim(name(1:idum-1))
                fraglib_sub(nfraglib)%s = trim(adjustl(name(idum+1:)))
             else
                fraglib_cat(nfraglib)%s = name
                fraglib_sub(nfraglib)%s = ""
             end if
             if (len_trim(fraglib_cat(nfraglib)%s) == 0) &
                fraglib_cat(nfraglib)%s = "Other"
          end if
       end if
    end do
    call fclose(lu)

    if (nfraglib > 0) then
       call realloc(fraglib_name,nfraglib)
       call realloc(fraglib_anchor,nfraglib)
       call realloc(fraglib_attach,nfraglib)
       call realloc(fraglib_cat,nfraglib)
       call realloc(fraglib_sub,nfraglib)
       call realloc(fraglib_radius,nfraglib)
    end if

  end subroutine fraglib_ensure

end submodule builder
