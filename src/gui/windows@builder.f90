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
  use types, only: vstring
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

  ! most neighbor directions considered when placing a ligand
  integer, parameter :: maxnbcand = 24

  ! bond styles in the local geometry pictograms
  integer, parameter :: sty_plain = 0 ! plain line, in the plane of the diagram
  integer, parameter :: sty_wedge = 1 ! filled wedge, towards the viewer
  integer, parameter :: sty_hash = 2  ! hashed wedge, away from the viewer

contains

  !> Draw the builder window. The builder is tied to its parent view
  !> (w%idparent) and operates on the system shown in it.
  module subroutine draw_builder(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_geometry,&
       reread_system_from_file, atlisttype_ncel_frac
    use dynamics, only: md_relax
    use energy, only: ff_backend_applicable, ff_backend_default
    use gui_main, only: ColorHighlightEditDistScene
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple, iw_dragfloat_real8,&
       iw_periodictable, iw_menuitem, invmult, mult
    use keybindings, only: is_bind_event, get_bind_keyname, BIND_PICKATOM_SELECT,&
       BIND_PICKATOM_ALT, BIND_RECALC_BONDS, BIND_NAV_MEASURE, BIND_EDIT_D_A_PHI,&
       BIND_REOPEN, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG, BIND_CANCEL
    use interfaces_glfw, only: glfwGetTime
    use tools_io, only: string, nameguess
    use tools_math, only: cross, axisangle2mat
    use param, only: bohrtoa, pi, eye
    class(window), intent(inout), target :: w

    integer :: isys, iview, icel, imode, nsel, iside, ibold, ibnew, izout, i, j, k
    logical :: doquit, goodparent, havesys, ok, lcommit, ldum, relaxing
    real*8 :: dist, ang
    real(c_float) :: xclick(2)
    character(len=:), allocatable :: stropt, errmsg
    character(kind=c_char,len=:), allocatable, target :: strcat, strsub

    logical, save :: ttshown = .false. ! tooltip flag

    real*8, parameter :: stale_gap = 1d0 ! discard clicks older than this (s)
    real*8, parameter :: eps_dzero = 1d-10 ! degenerate-distance threshold (bohr)
    ! bond order for each bond-combo index >= 1 (single, double, triple, dashed, aromatic)
    integer, parameter :: bondorder(5) = (/1,2,3,0,-1/)

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
       w%builder_vm = 0
       w%builder_isys = 0
       w%builder_time = 0d0
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

    ! handle an active builder mode commanded to the parent view
    if (w%builder_vm /= 0) then
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
             win(iview)%vmdata%tooltip_ig = w%builder_addatom_ig
             win(iview)%vmdata%tooltip_iz = w%builder_addatom_z
          elseif (allocated(w%builder_frag_name)) then
             win(iview)%vmdata%tooltip_frag = w%builder_frag_name
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
                if (w%builder_vm == vm_builder_addatom) then
                   call addatom_apply(xclick,icel)
                else
                   call addfrag_apply(xclick,icel)
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
          icel = win(iview)%vmdata%idx(1)
          imode = win(iview)%vmdata%flag
          win(iview)%vmdata%idx = 0
          if (glfwGetTime() - w%builder_time < stale_gap .and.&
             sysc(w%builder_isys)%timelastchange_geometry <= w%builder_time) then
             if (w%builder_vm == vm_builder_valence) then
                ! change valence: main pick = add a hydrogen, alternate
                ! pick = remove one
                call sysc(w%builder_isys)%change_valence(icel,imode)
             elseif (w%builder_vm == vm_builder_remove .and. imode == 1) then
                ! remove atoms (main pick only): the atom and its
                ! terminal hydrogens
                call sysc(w%builder_isys)%remove_atom_hydrogens(icel)
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

    ! header: the system this builder operates on
    call iw_text("System",highlight=.true.)
    if (havesys) &
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)
    if (iw_button("Restore",danger=.true.,disabled=.not.havesys)) then
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
       " first opened ("//trim(get_bind_keyname(BIND_REOPEN))//")",ttshown)

    ! the add section
    call iw_text("Add",highlight=.true.)
    if (w%builder_vm == vm_builder_addatom) then
       ! armed: the button cancels the mode
       if (iw_button("Add Atoms",danger=.true.)) &
          call builder_stop()
       call iw_tooltip("Stop adding atoms ("//trim(get_bind_keyname(BIND_CANCEL))//")",ttshown)
    else
       ldum = iw_button("Add Atoms",disabled=.not.havesys,popupcontext=ok,&
          popupflags=ImGuiPopupFlags_MouseButtonLeft)
       if (ok) then
          izout = iw_periodictable()
          if (izout > 0) then
             w%errmsg = ""
             w%builder_addatom_z = izout
             w%builder_addatom_ig = addatom_prefgeom(izout)
             call builder_toggle(vm_builder_addatom)
             call igCloseCurrentPopup()
          end if
          call igEndPopup()
       end if
       call iw_tooltip("Choose an element, then click in the view to add atoms of that"//&
          " element with the selected local geometry (the substituents are hydrogens)."//&
          " Clicking an atom replaces it with the new atom, keeping the substituents"//&
          " that are not terminal hydrogens",ttshown)
    end if
    ! the selected element, next to the button
    call iw_text(trim(nameguess(w%builder_addatom_z,.true.)),sameline=.true.)
    ! the local geometry selector: a grid of pictograms
    call iw_text("Local geometry",sameline=.true.)
    call draw_addatom_geom_grid(w%builder_addatom_ig)

    ! the add fragments section
    call iw_text("Add fragments",highlight=.true.)
    if (w%builder_vm == vm_builder_addfragment) then
       ! armed: the button cancels the mode
       if (iw_button("Add fragment",danger=.true.)) &
          call builder_stop()
       call iw_tooltip("Stop adding fragments ("//trim(get_bind_keyname(BIND_CANCEL))//")",ttshown)
    else
       ldum = iw_button("Add fragment",disabled=.not.havesys,popupcontext=ok,&
          popupflags=ImGuiPopupFlags_MouseButtonLeft)
       if (ok) then
          call fraglib_ensure()
          if (nfraglib == 0) then
             call iw_text("No fragments available")
          else
             ! one submenu per category, in order of first appearance, with
             ! the fragments of a "parent/child" category one level deeper
             do k = 1, nfraglib
                ! skip a category whose submenu has already been drawn
                if (.not.fraglib_isfirst(k,.false.)) cycle
                strcat = trim(fraglib_cat(k)%s) // c_null_char
                if (igBeginMenu(c_loc(strcat),.true._c_bool)) then
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
                      strsub = trim(fraglib_sub(j)%s) // c_null_char
                      if (igBeginMenu(c_loc(strsub),.true._c_bool)) then
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
          end if
          call igEndPopup()
       end if
       call iw_tooltip("Choose a fragment, then click in the view to add it.",ttshown)
    end if
    ! the selected fragment, next to the button
    if (allocated(w%builder_frag_name)) &
       call iw_text(trim(w%builder_frag_name),sameline=.true.)

    ! the atoms section
    call iw_text("Atoms",highlight=.true.)
    if (iw_button("Remove atoms",disabled=.not.havesys)) then
       w%errmsg = ""
       call builder_toggle(vm_builder_remove)
    end if
    call iw_tooltip("Remove ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//") the clicked"//&
       " atoms in the view, along with their terminal hydrogens",ttshown)

    ! the valence section
    call iw_text("Valence",highlight=.true.)
    if (iw_button("Change",disabled=.not.havesys)) then
       w%errmsg = ""
       call builder_toggle(vm_builder_valence)
    end if
    call iw_tooltip("Add ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//") or remove ("//&
       trim(get_bind_keyname(BIND_PICKATOM_ALT))//") hydrogens on the clicked atoms in the"//&
       " view, repositioning the terminal substituents",ttshown)

    ! the bonds section
    call iw_text("Bonds",highlight=.true.)
    if (iw_button("Rebond",disabled=.not.havesys)) then
       w%errmsg = ""
       call sysc(isys)%rebond()
    end if
    call iw_tooltip("Recompute the bond connectivity for this system ("//&
       trim(get_bind_keyname(BIND_RECALC_BONDS))//")",ttshown)

    ! number of measure-selected atoms in the parent view
    nsel = 0
    if (goodparent) then
       if (associated(win(iview)%sc)) nsel = win(iview)%sc%nmsel
    end if

    ! the keybinding request toggles an edit session: two selected
    ! atoms toggle edit distance, three edit angle, four edit dihedral,
    ! and otherwise deactivate the active session, if any
    if (w%edit_pending .or. (w%focused() .and. is_bind_event(BIND_EDIT_D_A_PHI))) then
       w%edit_pending = .false.
       if (nsel >= 2 .and. nsel <= 4) then
          call edit_toggle(nsel)
       elseif (w%edit_kind /= 0) then
          call w%edit_stop()
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

    ! edit distance section
    call iw_text("Edit distance",highlight=.true.)
    if (w%edit_kind /= 2) &
       call iw_text("Select two atoms in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//&
       "), then press Edit distance or "//trim(get_bind_keyname(BIND_EDIT_D_A_PHI)))
    ok = (w%edit_kind == 2) .or. (havesys .and. nsel == 2)
    if (iw_button("Edit distance",disabled=.not.ok)) then
       w%errmsg = ""
       call edit_toggle(2)
    end if
    call iw_tooltip("Edit the distance and bond between the two selected atoms ("//&
       trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//"); press again to stop",ttshown)

    if (w%edit_kind == 2) then
       ! the two atoms
       call edit_atom_labels()

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
       if (all(w%edit_imove == 0)) call igBeginDisabled(.true._c_bool)
       ldum = iw_dragfloat_real8("Distance (Å)##editdistdrag",x1=dist,speed=0.01d0,&
          min=0.1d0,scale=bohrtoa,decimal=3,notlive=.true.,committed=lcommit,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       if (all(w%edit_imove == 0)) then
          call igEndDisabled()
       else
          if (ldum) call editdist_apply(dist)
          if (lcommit) call edit_commit()
       end if
       call iw_tooltip("Distance between the two atoms (Å)",ttshown)

       ! apply button
       if (iw_button("Apply",danger=.true.)) then
          w%errmsg = ""
          call w%edit_stop()
       end if
       call iw_tooltip("Keep the current distance and release the two atoms",ttshown)
    end if

    ! edit angle section
    call iw_text("Edit angle",highlight=.true.)
    if (w%edit_kind /= 3) &
       call iw_text("Select three atoms in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//&
       "), then press Edit angle or "//trim(get_bind_keyname(BIND_EDIT_D_A_PHI)))
    ok = (w%edit_kind == 3) .or. (havesys .and. nsel == 3)
    if (iw_button("Edit angle",disabled=.not.ok)) then
       w%errmsg = ""
       call edit_toggle(3)
    end if
    call iw_tooltip("Edit the angle between the three selected atoms, vertex at the"//&
       " second atom ("//trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//"); press again to stop",ttshown)

    if (w%edit_kind == 3) then
       ! the three atoms
       call edit_atom_labels()

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
       if (all(w%edit_imove == 0)) call igBeginDisabled(.true._c_bool)
       ldum = iw_dragfloat_real8("Angle (°)##editangdrag",x1=ang,speed=0.2d0,&
          min=0d0,max=180d0,decimal=2,notlive=.true.,committed=lcommit,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       if (all(w%edit_imove == 0)) then
          call igEndDisabled()
       else
          if (ldum) call editang_apply(ang)
          if (lcommit) call edit_commit()
       end if
       call iw_tooltip("Angle between the three atoms (degrees)",ttshown)

       ! apply button
       if (iw_button("Apply##editang",danger=.true.)) then
          w%errmsg = ""
          call w%edit_stop()
       end if
       call iw_tooltip("Keep the current angle and release the three atoms",ttshown)
    end if

    ! edit dihedral section
    call iw_text("Edit dihedral",highlight=.true.)
    if (w%edit_kind /= 4) &
       call iw_text("Select four atoms in the view ("//trim(get_bind_keyname(BIND_NAV_MEASURE))//&
       "), then press Edit dihedral or "//trim(get_bind_keyname(BIND_EDIT_D_A_PHI)))
    ok = (w%edit_kind == 4) .or. (havesys .and. nsel == 4)
    if (iw_button("Edit dihedral",disabled=.not.ok)) then
       w%errmsg = ""
       call edit_toggle(4)
    end if
    call iw_tooltip("Edit the dihedral angle of the four selected atoms about the 2-3"//&
       " bond ("//trim(get_bind_keyname(BIND_EDIT_D_A_PHI))//"); press again to stop",ttshown)

    if (w%edit_kind == 4) then
       ! the four atoms
       call edit_atom_labels()

       ! terminal atom combos (atoms 2 and 3 always stay fixed)
       do iside = 1, 4, 3
          stropt = "Fixed"//c_null_char//"Rotate atom"//c_null_char
          if (w%edit_fragok(iside)) stropt = stropt // "Rotate group (fix 2-3)"//c_null_char
          ! fragok(2)/(3) = validity of the half adjacent to terminal 1/4
          if (w%edit_fragok(merge(2,3,iside == 1))) &
             stropt = stropt // "Rotate group (move 2-3)"//c_null_char
          call iw_combo_simple("Atom "//string(iside)//"##editdihmove"//string(iside),&
             stropt,w%edit_imove(iside))
          call iw_tooltip("What rotates about the 2-3 axis when the dihedral changes:"//&
             " nothing (fixed), the terminal atom, the group attached to it (the"//&
             " substituents of atoms 2 and 3 stay fixed), or the whole half severed"//&
             " at the 2-3 bond (the environments of atoms 2 and 3 stay rigid)",ttshown)
       end do

       ! dihedral drag-float (degrees)
       ang = editdih_val() * 180d0 / pi
       if (all(w%edit_imove == 0)) call igBeginDisabled(.true._c_bool)
       ldum = iw_dragfloat_real8("Dihedral (°)##editdihdrag",x1=ang,speed=0.2d0,&
          min=-180d0,max=180d0,decimal=2,notlive=.true.,committed=lcommit,&
          flags=ImGuiSliderFlags_AlwaysClamp)
       if (all(w%edit_imove == 0)) then
          call igEndDisabled()
       else
          if (ldum) call editdih_apply(ang)
          if (lcommit) call edit_commit()
       end if
       call iw_tooltip("Dihedral angle of the four atoms (degrees)",ttshown)

       ! apply button
       if (iw_button("Apply##editdih",danger=.true.)) then
          w%errmsg = ""
          call w%edit_stop()
       end if
       call iw_tooltip("Keep the current dihedral and release the four atoms",ttshown)
    end if

    ! transient highlight of the latched atoms
    if (w%edit_kind /= 0) &
       call sysc(w%edit_isys)%highlight_atoms(.true.,w%edit_idx(1,1:w%edit_kind),&
       atlisttype_ncel_frac,spread(ColorHighlightEditDistScene,2,w%edit_kind))

    ! the relax section
    call iw_text("Relax",highlight=.true.)
    relaxing = .false.
    if (havesys) then
       relaxing = sysc(isys)%md_run .and. sysc(isys)%md%mode == md_relax
       ! resolve the method before the button (the combo below would
       ! only snap an inapplicable backend after this frame's click)
       if (sysc(isys)%md_backend < 0 .or.&
          .not.ff_backend_applicable(sysc(isys)%md_backend,sys(isys)%c)) &
          sysc(isys)%md_backend = ff_backend_default(sys(isys)%c)
    end if
    if (.not.relaxing) then
       ok = havesys
       if (ok) ok = .not.sysc(isys)%md_run
       if (iw_button("Relax",danger=.true.,disabled=.not.ok)) then
          ! release any edit session (a running relaxation ends it anyway)
          call w%edit_stop()
          call sysc(isys)%md_start(md_relax,errmsg)
          w%errmsg = errmsg
          if (len_trim(w%errmsg) == 0) then
             if (associated(win(iview)%sc)) win(iview)%forcerender = .true.
          end if
       end if
       call iw_tooltip("Relax the geometry to the nearest energy minimum, stopping when"//&
          " the maximum force falls below the threshold",ttshown)
    else
       if (iw_button("Stop",danger=.true.)) then
          call sysc(isys)%md_stop()
          if (associated(win(iview)%sc)) win(iview)%forcerender = .true.
       end if
       call iw_tooltip("Stop the geometry relaxation ("//&
          trim(get_bind_keyname(BIND_CANCEL))//")",ttshown)
    end if

    if (havesys) then
       ! method combo, next to the button
       call draw_ff_backend_combo(isys,"##builderrelaxmethod",15)
       call iw_tooltip("Method for the calculation of energies and forces",ttshown)

       ! force convergence threshold
       ldum = iw_dragfloat_real8("Force threshold (eV/Å)##relaxfconv",x1=sysc(isys)%md%fconv,&
          speed=0.0005d0,min=0.0001d0,max=1d0,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Stop the relaxation when the maximum force falls below this value",ttshown)

       ! live status: energy and maximum force
       if (sysc(isys)%md%ready .and. sysc(isys)%md%nat > 0 .and.&
          sysc(isys)%md%mode == md_relax) then
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

    ! the symmetry section
    call iw_text("Symmetry",highlight=.true.)
    if (iw_button("Refine symmetry",disabled=.not.havesys)) then
       w%errmsg = ""
       call sysc(isys)%refine_symmetry(w%errmsg)
       ! hold the camera on the parent view's scene (the main scene when
       ! the parent is the main view, its private scene otherwise)
       if (associated(win(iview)%sc)) win(iview)%sc%nextbuildlists_fixcam = .true.
    end if
    call iw_tooltip("Refine the geometry having atoms occupy their symmetry positions exactly",ttshown)

    ! error message, if any
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

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

    ! Show the labels of the latched atoms in one row.
    subroutine edit_atom_labels()
      integer :: is

      do is = 1, w%edit_kind
         call iw_text("("//string(is)//") "//&
            trim(sys(isys)%c%at(sys(isys)%c%atcel(w%edit_idx(1,is))%idx)%name)//&
            string(w%edit_idx(1,is)),sameline=(is > 1))
      end do
    end subroutine edit_atom_labels

    ! Rotate the moving terminal atoms (or their groups) about the
    ! axis through center: the first latched atom rotates by
    ! -dang/nrot and the last by +dang/nrot, which opens the angle or
    ! dihedral by dang. In the dihedral "move 2-3" mode (imove = 3)
    ! the move set is the corresponding half, stored in the interior
    ! columns.
    subroutine edit_rotate_terminals(rnew,center,axis,dang,nrot)
      real*8, intent(inout) :: rnew(:,:)
      real*8, intent(in) :: center(3), axis(3), dang
      integer, intent(in) :: nrot

      integer :: i, k, n, icol
      real*8 :: rot(3,3)

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
            rnew(:,w%edit_frag(k,icol)) = center + matmul(rot,rnew(:,w%edit_frag(k,icol)) - center)
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
      integer, allocatable :: ifr(:), work(:,:)

      allocate(work(sys(isys)%c%ncel,ikind))
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
               w%edit_idx(1,3),n,ifr,disc)
         elseif (is > 1 .and. is < ikind) then
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is-1),n,ifr,disc,w%edit_idx(1,is),w%edit_idx(1,is+1))
         elseif (is > 1) then
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is-1),n,ifr,disc)
         else
            it = is
            call sys(isys)%c%masked_fragment(w%edit_idx(1,is),w%edit_idx(1,is),&
               w%edit_idx(1,is+1),n,ifr,disc)
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
      end do
      if (ikind == 4) then
         ! offer the "move 2-3" mode only when the plain group mode is
         ! also available: the combo option lists are built contiguously,
         ! so entry 3 (move) must not exist without entry 2 (fix). This
         ! can suppress a valid half when the terminal's own group is
         ! invalid (e.g. a ring through the terminal-axis bond).
         if (.not.w%edit_fragok(1)) w%edit_fragok(2) = .false.
         if (.not.w%edit_fragok(4)) w%edit_fragok(3) = .false.
      end if
      if (allocated(w%edit_frag)) deallocate(w%edit_frag)
      allocate(w%edit_frag(maxval(w%edit_nfrag(1:ikind)),ikind))
      w%edit_frag = 0
      do is = 1, ikind
         w%edit_frag(1:w%edit_nfrag(is),is) = work(1:w%edit_nfrag(is),is)
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
      integer :: ia, j, k, lv(3)

      idx = 0
      ia = w%edit_idx(1,1)
      lv = w%edit_idx(2:4,2) - w%edit_idx(2:4,1)
      if (.not.allocated(sys(isys)%c%nstar)) return
      do j = 1, sys(isys)%c%nstar(ia)%ncon
         if (sys(isys)%c%nstar(ia)%idcon(j) == w%edit_idx(1,2) .and.&
            all(sys(isys)%c%nstar(ia)%lcon(:,j) == lv)) then
            do k = 1, size(bondorder,1)
               if (sys(isys)%c%nstar(ia)%ordcon(j) == bondorder(k)) then
                  idx = k
                  return
               end if
            end do
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

    ! Unit vector perpendicular to v (assumed nonzero).
    function editang_perp(v) result(u)
      real*8, intent(in) :: v(3)
      real*8 :: u(3)

      if (abs(v(1)) < 0.9d0 * norm2(v)) then
         u = cross(v,(/1d0,0d0,0d0/))
      else
         u = cross(v,(/0d0,1d0,0d0/))
      end if
      u = u / norm2(u)
    end function editang_perp

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
            w0 = editang_perp(va)
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
         if (norm2(axis) < eps_dzero) axis = editang_perp(va)
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

      if (goodparent) call win(iview)%viewmode_release_forced(w%id,w%builder_vm)
      w%builder_vm = 0
      w%builder_isys = 0
    end subroutine builder_stop

    ! Unproject the clicked texture position xpos onto the plane
    ! parallel to the screen through the scene center: x0 is the
    ! clicked point in world coordinates and the columns of rcam are
    ! the camera axes (screen right, screen up, towards the viewer).
    subroutine builder_click_frame(xpos,x0,rcam,ok)
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

    end subroutine builder_click_frame

    ! Add-atoms mode: on empty space (icel = 0), unproject the clicked
    ! texture position xpos to the plane parallel to the screen
    ! through the scene center and add the fragment there. On an atom
    ! (icel > 0), replace it (see addatom_replace).
    subroutine addatom_apply(xpos,icel)
      real(c_float), intent(in) :: xpos(2)
      integer, intent(in) :: icel

      real*8 :: x0(3), rot(3,3), bl
      real*8 :: tvec(3,maxaddsub), xat(3,maxaddsub+1)
      integer :: zat(maxaddsub+1), nsub, nat, k
      logical :: ok

      call builder_click_frame(xpos,x0,rot,ok)
      if (.not.ok) return

      ! clicked on an atom: replace it instead
      if (icel > 0) then
         call addatom_replace(icel,rot)
         return
      end if

      ! the fragment: central atom plus hydrogens (per-system covalent
      ! radii, same source as change_valence)
      call addatom_template(w%builder_addatom_ig,nsub,tvec)
      nat = 1 + nsub
      zat(1) = w%builder_addatom_z
      xat(:,1) = x0
      if (nsub > 0) then
         bl = sysc(w%builder_isys)%atmcov(zat(1)) + sysc(w%builder_isys)%atmcov(1)
         do k = 1, nsub
            zat(1+k) = 1
            xat(:,1+k) = x0 + bl * matmul(rot,tvec(:,k))
         end do
      end if
      call sysc(w%builder_isys)%add_atoms_fragment(nat,zat(1:nat),xat(:,1:nat))

    end subroutine addatom_apply

    ! Add-atoms mode, click on cell atom icel: replace it with the fragment.
    subroutine addatom_replace(icel,rcam)
      integer, intent(in) :: icel
      real*8, intent(in) :: rcam(3,3)

      integer :: isysl, nsub, m, ndel, nadd, k, nb, ncon, znew, zanchor
      integer :: zat(maxaddsub+1)
      integer, allocatable :: idel(:)
      real*8 :: tvec(3,maxaddsub), xat(3,maxaddsub+1), rot(3,3)
      real*8 :: x0(3), dir(3), bl, dnorm
      real*8, allocatable :: ufix(:,:)
      logical :: used(maxaddsub), isterm

      isysl = w%builder_isys
      if (icel < 1 .or. icel > sys(isysl)%c%ncel) return
      if (.not.allocated(sys(isysl)%c%nstar)) return
      znew = w%builder_addatom_z
      call addatom_template(w%builder_addatom_ig,nsub,tvec)

      ! classify the neighbors: kept substituents vs terminal hydrogens
      x0 = sys(isysl)%c%atcel(icel)%r
      ncon = sys(isysl)%c%nstar(icel)%ncon
      allocate(idel(ncon+1),ufix(3,max(ncon,1)))
      ndel = 1
      idel(1) = icel
      m = 0
      zanchor = 0
      do k = 1, ncon
         nb = sys(isysl)%c%nstar(icel)%idcon(k)
         isterm = ncon > 1 .and. sys(isysl)%c%spc(sys(isysl)%c%atcel(nb)%is)%z == 1
         if (isterm) isterm = sys(isysl)%c%nstar(nb)%ncon == 1
         if (isterm) then
            if (.not.any(idel(1:ndel) == nb)) then
               ndel = ndel + 1
               idel(ndel) = nb
            end if
         else
            m = m + 1
            ufix(:,m) = sys(isysl)%c%x2c(sys(isysl)%c%atcel(nb)%x +&
               dble(sys(isysl)%c%nstar(icel)%lcon(:,k)))
            zanchor = sys(isysl)%c%spc(sys(isysl)%c%atcel(nb)%is)%z
         end if
      end do

      ! on a single kept substituent, re-bond the new atom at the
      ! covalent distance along the old bond direction
      if (m == 1) then
         dir = x0 - ufix(:,1)
         dnorm = norm2(dir)
         if (dnorm < eps_dzero) return
         x0 = ufix(:,1) + (sysc(isysl)%atmcov(zanchor) + sysc(isysl)%atmcov(znew)) *&
            dir / dnorm
      end if

      ! enough kept substituents already: replace the element and remove the terminal hydrogens
      if (m >= nsub) then
         zat(1) = znew
         xat(:,1) = x0
         call sysc(isysl)%replace_atoms_fragment(ndel,idel(1:ndel),1,zat(1:1),xat(:,1:1))
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
      bl = sysc(isysl)%atmcov(znew) + sysc(isysl)%atmcov(1)
      do k = 1, nsub
         if (used(k)) cycle
         nadd = nadd + 1
         zat(nadd) = 1
         xat(:,nadd) = x0 + bl * matmul(rot,tvec(:,k))
      end do
      call sysc(isysl)%replace_atoms_fragment(ndel,idel(1:ndel),nadd,zat(1:nadd),xat(:,1:nadd))

    end subroutine addatom_replace

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
      integer, intent(in) :: i

      if (iw_menuitem(trim(fraglib_name(i)%s))) then
         w%errmsg = ""
         if (addfrag_load(i)) call builder_toggle(vm_builder_addfragment)
         call igCloseCurrentPopup()
      end if

    end subroutine fragment_menuitem

    ! Load fragment ifrag of the library into the builder window state.
    ! Returns .false. (and sets w%errmsg) if the fragment cannot be read
    ! or its anchor/attach atoms are not valid.
    function addfrag_load(ifrag) result(ok)
      use crystalseedmod, only: crystalseed
      use global, only: flib_file
      use tools_io, only: string
      integer, intent(in) :: ifrag
      logical :: ok

      type(crystalseed) :: seed
      integer :: i, ia, it

      ok = .false.
      w%builder_frag_nat = 0
      if (ifrag < 1 .or. ifrag > nfraglib) return

      ! read the coordinates with the library reader
      call seed%read_library(trim(fraglib_name(ifrag)%s),ok,file=flib_file)
      if (.not.ok .or. seed%nat <= 0 .or. seed%nspc <= 0) then
         w%errmsg = "Could not read fragment " // trim(fraglib_name(ifrag)%s) //&
            " from the fragment library"
         ok = .false.
         return
      end if

      ! the anchor and its attachment placeholder must exist
      ia = fraglib_anchor(ifrag)
      it = fraglib_attach(ifrag)
      ok = (ia >= 1 .and. ia <= seed%nat .and. it >= 1 .and. it <= seed%nat .and. ia /= it)
      if (.not.ok) then
         w%errmsg = "Fragment " // trim(fraglib_name(ifrag)%s) //&
            " has an invalid anchor or attach atom"
         return
      end if

      ! save the fragment (seed%x is Cartesian and in bohr for a molecule)
      if (allocated(w%builder_frag_z)) deallocate(w%builder_frag_z)
      if (allocated(w%builder_frag_x)) deallocate(w%builder_frag_x)
      allocate(w%builder_frag_z(seed%nat),w%builder_frag_x(3,seed%nat))
      do i = 1, seed%nat
         w%builder_frag_z(i) = seed%spc(seed%is(i))%z
         w%builder_frag_x(:,i) = seed%x(:,i)
      end do
      w%builder_frag_nat = seed%nat
      w%builder_frag_ianchor = ia
      w%builder_frag_iattach = it
      w%builder_frag_radius = fraglib_radius(ifrag)
      w%builder_frag_name = trim(fraglib_name(ifrag)%s)

    end function addfrag_load

    ! Add-fragments mode: place the loaded fragment. On empty space
    ! (icel = 0) the fragment goes at the clicked position, camera-aligned
    ! and capped with its placeholder atom. On an atom (icel > 0) the
    ! anchor replaces it (see addfrag_target), unless the fragment is a
    ! ligand, in which case it bonds to it (see addfrag_ligand_target).
    subroutine addfrag_apply(xpos,icel)
      integer, intent(in) :: icel
      real(c_float), intent(in) :: xpos(2)

      integer :: isysl, i, ia, nadd, ndel
      integer, allocatable :: idel(:), zadd(:)
      real*8, allocatable :: xadd(:,:)
      real*8 :: x0(3), rcam(3,3), e(3,3), rot(3,3), p0(3), t1(3), dattach(3)
      logical :: ok, keepattach

      isysl = w%builder_isys
      ia = w%builder_frag_ianchor
      if (w%builder_frag_nat <= 0) return
      call builder_click_frame(xpos,x0,rcam,ok)
      if (.not.ok) return

      ! the fragment frame and where the fragment goes
      call addfrag_frame(e)
      ! the placeholder caps the anchor on empty space, but only if it is
      ! a real atom: a dummy (the dative ligands) marks a direction only
      keepattach = (icel == 0 .and. w%builder_frag_z(w%builder_frag_iattach) > 0)
      ndel = 0
      allocate(idel(1))
      if (icel == 0) then
         ! empty space: the attachment points to the left of the screen
         p0 = x0
         dattach = -rcam(:,1)
      elseif (w%builder_frag_radius > 0d0) then
         ! a ligand: it bonds to the clicked atom instead of replacing it
         call addfrag_ligand_target(icel,rcam,p0,dattach,ok)
         if (.not.ok) return
      else
         call addfrag_target(icel,rcam,p0,dattach,ndel,idel,ok)
         if (.not.ok) return
      end if

      ! orient the fragment: the attachment along dattach and the plane
      ! of the fragment facing the viewer as much as it allows
      call addfrag_rotation(e,dattach,real(rcam(:,3),8),rot)
      t1 = matmul(rot,e(:,1))

      ! the fragment atoms, rotated and translated onto the anchor; the
      ! placeholder is dropped when the anchor bonds to the structure
      allocate(zadd(w%builder_frag_nat),xadd(3,w%builder_frag_nat))
      nadd = 0
      do i = 1, w%builder_frag_nat
         if (.not.keepattach .and. i == w%builder_frag_iattach) cycle
         nadd = nadd + 1
         zadd(nadd) = w%builder_frag_z(i)
         xadd(:,nadd) = p0 + matmul(rot,w%builder_frag_x(:,i) - w%builder_frag_x(:,ia))
      end do

      ! cap the anchor at the sum of covalent radii
      if (keepattach) &
         xadd(:,w%builder_frag_iattach) = p0 + t1 * &
         (sysc(isysl)%atmcov(w%builder_frag_z(ia)) +&
         sysc(isysl)%atmcov(w%builder_frag_z(w%builder_frag_iattach)))

      call sysc(isysl)%replace_atoms_fragment(ndel,idel(1:max(ndel,1)),nadd,&
         zadd(1:nadd),xadd(:,1:nadd))

    end subroutine addfrag_apply

    ! Add-fragments mode, click on cell atom icel: the anchor of the
    ! fragment replaces it. Returns the position of the anchor (p0), the
    ! direction from the anchor towards the atoms it bonds to (dattach),
    ! and the atoms to delete (icel plus the terminal hydrogens that make
    ! room for the fragment).
    subroutine addfrag_target(icel,rcam,p0,dattach,ndel,idel,ok)
      integer, intent(in) :: icel
      real*8, intent(in) :: rcam(3,3)
      real*8, intent(out) :: p0(3), dattach(3)
      integer, intent(out) :: ndel
      integer, allocatable, intent(inout) :: idel(:)
      logical, intent(out) :: ok

      integer :: isysl, ncon, k, nb, nanch, nrem, m, zkept
      real*8 :: x0(3), xnb(3), xkept(3), dir(3), dn, sumdir(3)
      logical :: isterm

      ok = .false.
      isysl = w%builder_isys
      if (icel < 1 .or. icel > sys(isysl)%c%ncel) return
      if (.not.allocated(sys(isysl)%c%nstar)) return

      x0 = sys(isysl)%c%atcel(icel)%r
      ncon = sys(isysl)%c%nstar(icel)%ncon
      if (allocated(idel)) deallocate(idel)
      allocate(idel(ncon+1))
      ndel = 1
      idel(1) = icel

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
      xkept = 0d0
      do k = 1, ncon
         nb = sys(isysl)%c%nstar(icel)%idcon(k)
         isterm = (sys(isysl)%c%spc(sys(isysl)%c%atcel(nb)%is)%z == 1)
         if (isterm) isterm = (sys(isysl)%c%nstar(nb)%ncon == 1)
         if (isterm .and. nrem < nanch) then
            if (any(idel(1:ndel) == nb)) cycle
            nrem = nrem + 1
            ndel = ndel + 1
            idel(ndel) = nb
         else
            xnb = sys(isysl)%c%x2c(sys(isysl)%c%atcel(nb)%x +&
               dble(sys(isysl)%c%nstar(icel)%lcon(:,k)))
            dir = xnb - x0
            dn = norm2(dir)
            if (dn < eps_dzero) cycle
            m = m + 1
            sumdir = sumdir + dir / dn
            xkept = xnb
            zkept = sys(isysl)%c%spc(sys(isysl)%c%atcel(nb)%is)%z
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
            dir = x0 - xkept
            dn = norm2(dir)
            if (dn < eps_dzero) return
            p0 = xkept + (sysc(isysl)%atmcov(zkept) +&
               sysc(isysl)%atmcov(w%builder_frag_z(w%builder_frag_ianchor))) * dir / dn
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

      integer :: isysl, ncon, k, nb, m, i1, i2, i3, ncand, ibest
      real*8 :: x0(3), xnb(3), dir(3), dn, sumdir(3), smax, sbest
      real*8 :: u(3,maxnbcand), cand(3,27)

      ok = .false.
      isysl = w%builder_isys
      if (icel < 1 .or. icel > sys(isysl)%c%ncel) return
      if (.not.allocated(sys(isysl)%c%nstar)) return

      ! unit directions to the neighbors of the clicked atom
      x0 = sys(isysl)%c%atcel(icel)%r
      ncon = sys(isysl)%c%nstar(icel)%ncon
      sumdir = 0d0
      m = 0
      do k = 1, ncon
         if (m >= maxnbcand) exit
         nb = sys(isysl)%c%nstar(icel)%idcon(k)
         xnb = sys(isysl)%c%x2c(sys(isysl)%c%atcel(nb)%x +&
            dble(sys(isysl)%c%nstar(icel)%lcon(:,k)))
         dir = xnb - x0
         dn = norm2(dir)
         if (dn < eps_dzero) cycle
         m = m + 1
         u(:,m) = dir / dn
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
         w%builder_frag_radius) * dattach
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

      ia = w%builder_frag_ianchor
      n = 1
      do i = 1, w%builder_frag_nat
         if (i == ia .or. i == w%builder_frag_iattach) cycle
         if (norm2(w%builder_frag_x(:,i) - w%builder_frag_x(:,ia)) < bondfactor *&
            (sysc(w%builder_isys)%atmcov(w%builder_frag_z(i)) +&
            sysc(w%builder_isys)%atmcov(w%builder_frag_z(ia)))) n = n + 1
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

      ! e1: the attachment direction
      e(:,1) = w%builder_frag_x(:,w%builder_frag_iattach) -&
         w%builder_frag_x(:,w%builder_frag_ianchor)
      e(:,1) = e(:,1) / norm2(e(:,1))

      ! e3: the direction in which the fragment is least spread out
      cov = 0d0
      do i = 1, w%builder_frag_nat
         xd = w%builder_frag_x(:,i) - w%builder_frag_x(:,w%builder_frag_ianchor)
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
         e(:,3) = (/1d0,0d0,0d0/) - e(1,1) * e(:,1)
         dn = norm2(e(:,3))
         if (dn < eps_dzero) then
            e(:,3) = (/0d0,1d0,0d0/) - e(2,1) * e(:,1)
            dn = norm2(e(:,3))
         end if
      end if
      e(:,3) = e(:,3) / dn
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
         t(:,3) = (/1d0,0d0,0d0/) - t(1,1) * t(:,1)
         dn = norm2(t(:,3))
         if (dn < eps_dzero) then
            t(:,3) = (/0d0,1d0,0d0/) - t(2,1) * t(:,1)
            dn = norm2(t(:,3))
         end if
      end if
      t(:,3) = t(:,3) / dn
      t(:,2) = cross(t(:,3),t(:,1))
      rot = matmul(t,transpose(e))

    end subroutine addfrag_rotation



    ! Toggle builder mode jvm (a vm_* constant): start it on the parent
    ! view, or stop it if it is already the active mode. Starting a mode
    ! while another one is active switches to the new mode.
    subroutine builder_toggle(jvm)
      integer, intent(in) :: jvm

      character(len=:), allocatable :: msg

      if (w%builder_vm == jvm) then
         call builder_stop()
      else
         if (w%builder_vm /= 0) call builder_stop()
         w%builder_vm = jvm
         w%builder_isys = isys
         w%builder_time = glfwGetTime()
         if (jvm == vm_builder_valence) then
            msg = "Add ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//&
               ") or remove ("//trim(get_bind_keyname(BIND_PICKATOM_ALT))//&
               ") hydrogens..."
         elseif (jvm == vm_builder_addatom) then
            msg = "Add atoms ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//&
               ") at the clicked positions or replace the clicked atoms"
         elseif (jvm == vm_builder_addfragment) then
            msg = "Add "//trim(w%builder_frag_name)//" ("//&
               trim(get_bind_keyname(BIND_PICKATOM_SELECT))//&
               ") at the clicked positions or "
            if (w%builder_frag_radius > 0d0) then
               msg = msg // "bond it to the clicked atom"
            else
               msg = msg // "replace the clicked atoms"
            end if
         else
            msg = "Remove ("//trim(get_bind_keyname(BIND_PICKATOM_SELECT))//&
               ") the atoms and their hydrogens"
         end if
         call win(iview)%viewmode_set_forced(jvm,msg,w%id)
      end if
    end subroutine builder_toggle
  end subroutine draw_builder

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
  !> index), updated when an icon is clicked. Hovering an icon shows
  !> the geometry name.
  subroutine draw_addatom_geom_grid(ig)
    use gui_main, only: g
    use utils, only: iw_tooltip
    integer, intent(inout) :: ig

    integer :: i
    logical :: hovered
    type(ImVec2) :: sz, p0, p1
    type(c_ptr) :: dl
    integer(c_int) :: col
    real(c_float) :: side
    character(len=:,kind=c_char), allocatable, target :: str1

    integer, parameter :: ncol = 7 ! icons per row

    side = 2.2_c_float * igGetTextLineHeightWithSpacing()
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
    real(c_float) :: ux, uy, blen, f, hw, rstart, sluma
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

    do k = 1, nb
       ! bond direction on screen (angles counterclockwise from +x, screen y points down)
       a = ang(k) * pi / 180d0
       ux = real(cos(a),c_float)
       uy = -real(sin(a),c_float)
       blen = rstart + (blrat * real(lfac(k),c_float) - r0rat) * side
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
