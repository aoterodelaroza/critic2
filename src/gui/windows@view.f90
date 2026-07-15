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

! Windows, view.
submodule (windows) view
  use interfaces_cimgui
  implicit none

  ! 4-identity matrix (c-float)
  real(c_float), parameter :: zero = 0._c_float
  real(c_float), parameter :: one = 1._c_float
  real(c_float), parameter :: eye4(4,4) = reshape((/&
     one,zero,zero,zero,&
     zero,one,zero,zero,&
     zero,zero,one,zero,&
     zero,zero,zero,one/),shape(eye4))

  ! ilock parameters for mouse interaction with view
  integer, parameter :: ilock_no = 0
  integer, parameter :: ilock_left = 1
  integer, parameter :: ilock_right = 2
  integer, parameter :: ilock_scroll = 3
  integer, parameter :: ilock_middle = 4
  integer, parameter :: ilock_mddrag = 5    ! dragging an atom during interactive dynamics
  integer, parameter :: ilock_mdmovemol = 6 ! rigidly translating a molecule during interactive dynamics
  integer, parameter :: ilock_mdrotmol = 7  ! rigidly rotating a molecule during interactive dynamics

  ! minimum time elapsed between consecutive queries of the pick buffer (seconds)
  real*8, parameter :: pick_interval = 1d0 / 10d0

  ! render-texture sizing (pixels). The texture is sized up to the view window
  ! (snapped to a bucket, capped) so the still image stays sharp; while the user
  ! manipulates the camera the scene is rendered at interactive_texture_side to
  ! keep interaction fast regardless of window size.
  integer(c_int), parameter :: min_texture_side = 1024_c_int
  integer(c_int), parameter :: max_texture_side = 4096_c_int
  integer(c_int), parameter :: interactive_texture_side = 1024_c_int
  integer(c_int), parameter :: texture_side_bucket = 256_c_int

contains

  !> Draw the view.
  module subroutine draw_view(w)
    use interfaces_glfw, only: glfwGetTime
    use interfaces_opengl3
    use interfaces_cimgui
    use keybindings, only: is_bind_event, BIND_VIEW_INC_NCELL, BIND_VIEW_DEC_NCELL,&
       BIND_VIEW_ALIGN_A_AXIS, BIND_VIEW_ALIGN_B_AXIS, BIND_VIEW_ALIGN_C_AXIS,&
       BIND_VIEW_ALIGN_X_AXIS, BIND_VIEW_ALIGN_Y_AXIS, BIND_VIEW_ALIGN_Z_AXIS,&
       BIND_VIEW_TOGGLE_ATOMS, BIND_VIEW_TOGGLE_BONDS, BIND_VIEW_CYCLE_LABELS,&
       BIND_VIEW_TOGGLE_CELL, BIND_VIEW_TOGGLE_POLYHEDRA, BIND_RECALC_BONDS,&
       get_bind_keyname, BIND_EDITSELECT_REMOVE, BIND_EDITSELECT_DESELECT,&
       BIND_EDITSELECT_SELECT_ALL, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use representations, only: reptype_atoms, reptype_unitcell, reptype_axes, reptype_symelem,&
       repflavor_atoms_ballandstick, repflavor_atoms_criticalpoints, repflavor_atoms_gradientpaths,&
       repflavor_atoms_vdwcontacts, repflavor_atoms_hbonds,&
       repflavor_atoms_sticks, repflavor_atoms_licorice, repflavor_unitcell_basic,&
       repflavor_axes, repflavor_atoms_polyhedra, repflavor_symelem
    use utils, only: iw_calcheight, iw_calcwidth, iw_clamp_color3, iw_combo_simple,&
       iw_setposx_fromend, iw_checkbox, iw_coloredit, iw_menuitem, iw_dragfloat_realc,&
       iw_text, iw_button, iw_tooltip, iw_intstepper, iw_radiobutton
    use crystalmod, only: iperiod_vacthr
    use global, only: dunit0, iunit_ang
    use systems, only: sysc, sys, sys_init, nsys, ok_system, are_threads_running
    use gui_main, only: g, fontsize, lockbehavior, tree_select_updates_view,&
       ColorBlack, ColorWhite, ColorClearTransparent
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: i, j, k, nrep, is, icel, ineq, iaux
    type(ImVec2) :: szavail, sz0, sz1, szero, pos
    type(ImVec4) :: tintcol, bgcol
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    character(len=:), allocatable, target :: msg
    logical(c_bool) :: is_selected
    logical :: hover, chbuild, chrender, goodsys, ldum, ok, ismol, isatom, isbond
    logical :: isuc, islabelsl, needpick, enabled, ispoly
    integer :: islabels
    logical :: ch
    integer(c_int) :: flags, nc(3), ires, idum
    integer(c_int) :: newside, vside
    real(c_float) :: scal, width, rgba(4)
    real(c_float) :: rscale, tmpuv
    logical :: interacting, selcleared
    real*8 :: x0(3), time
    type(ImVec2) :: sz
    logical :: changedisplay(5) ! 1=atoms, 2=bonds, 3=labels, 4=cell, 5=polyhedra

    logical, save :: ttshown = .false. ! tooltip flag

    ! coordinate this with objects (representation) menu in scenes module
    integer(c_int), parameter :: ic_closebutton = 0
    integer(c_int), parameter :: ic_viewbutton = 1
    integer(c_int), parameter :: ic_name = 2
    integer(c_int), parameter :: ic_type = 3
    integer(c_int), parameter :: ic_editbutton = 4

    ! initialize
    time = glfwGetTime()
    szero%x = 0
    szero%y = 0
    chrender = .false.
    chbuild = .false.
    if (w%firstpass) then
       w%mousepos_lastpick%x = 0._c_float
       w%mousepos_lastpick%y = 0._c_float
       w%viewmode = vm_navigate
       w%viewmode_transient = .false.
    end if

    !! update the tree based on time signals between dependent windows
    ! track the tree system if this is the main view and the option is active
    if (w%timelast_view_assign < win(iwin_tree)%timelast_tree_assign .and.&
       w%ismain .and. tree_select_updates_view) then
       call w%select_view(win(iwin_tree)%tree_selected)
    end if

    ! whether the selected view system is a good system, and associate the scene
    goodsys = ok_system(w%view_selected,sys_init)
    if (goodsys) then
       if (w%ismain) then
          if (.not.associated(w%sc)) w%sc => sysc(w%view_selected)%sc
       else
          if (.not.associated(w%sc)) allocate(w%sc)
          if (w%sc%isinit == 0) call w%sc%init(w%view_selected)
       end if
    else
       if (w%ismain) then
          nullify(w%sc)
       else
          if (associated(w%sc)) then
             call w%sc%end()
             deallocate(w%sc)
          end if
          w%forcerender = .true.
       end if
    end if

    ! flags for shortcuts
    isatom = .false.
    isbond = .false.
    islabels = -1
    islabelsl = .false.
    isuc = .false.
    ispoly = .false.
    if (associated(w%sc)) then
       do i = 1, w%sc%nrep
          if (w%sc%rep(i)%isinit) then
             if (w%sc%rep(i)%type == reptype_atoms) then
                isatom = isatom .or. w%sc%rep(i)%atoms%display
                isbond = isbond .or. w%sc%rep(i)%bonds%display
                ispoly = ispoly .or. w%sc%rep(i)%poly%display
                if (w%sc%rep(i)%labels%display .and. w%sc%rep(i)%flavor/=repflavor_atoms_criticalpoints .and.&
                   w%sc%rep(i)%flavor/=repflavor_atoms_gradientpaths) &
                   islabels = w%sc%rep(i)%labels%type
             elseif (w%sc%rep(i)%type == reptype_unitcell) then
                isuc = isuc .or. w%sc%rep(i)%shown
             end if
          end if
       end do
       if (islabels >= 0) islabelsl = .true.
    end if

    ! process shortcut display bindings
    if (associated(w%sc)) then
       changedisplay = .false.
       if (is_bind_event(BIND_VIEW_TOGGLE_ATOMS)) then
          changedisplay(1) = .true.
          isatom = .not.isatom
       end if
       if (is_bind_event(BIND_VIEW_TOGGLE_BONDS)) then
          changedisplay(2) = .true.
          isbond = .not.isbond
       end if
       if (is_bind_event(BIND_VIEW_CYCLE_LABELS)) then
          changedisplay(3) = .true.
          if (islabels == -1) then ! no labels -> atom-name
             islabels = 1
          elseif (islabels == 1) then ! atom-name -> celatom
             islabels = 2
          elseif (islabels == 2 .and..not.sys(w%view_selected)%c%ismolecule) then ! celatom -> wyckoff
             islabels = 8
          else ! * -> no labels
             islabels = -1
          end if
          islabelsl = (islabels >= 0)
       end if
       if (is_bind_event(BIND_VIEW_TOGGLE_CELL)) then
          changedisplay(4) = .true.
          isuc = .not.isuc
       end if
       if (is_bind_event(BIND_VIEW_TOGGLE_POLYHEDRA)) then
          changedisplay(5) = .true.
          ispoly = .not.ispoly
       end if
       if (is_bind_event(BIND_RECALC_BONDS)) &
          call sysc(w%view_selected)%rebond()
       if (any(changedisplay)) then
          do i = 1, w%sc%nrep
             if (w%sc%rep(i)%isinit) then
                if (w%sc%rep(i)%type == reptype_atoms) then
                   if (changedisplay(1)) w%sc%rep(i)%atoms%display = isatom
                   if (changedisplay(2)) w%sc%rep(i)%bonds%display = isbond
                   if (changedisplay(5)) w%sc%rep(i)%poly%display = ispoly
                   if (changedisplay(3) .and. w%sc%rep(i)%flavor/=repflavor_atoms_criticalpoints .and.&
                      w%sc%rep(i)%flavor/=repflavor_atoms_gradientpaths) then
                      w%sc%rep(i)%labels%display = islabelsl
                      if (islabelsl) then
                         w%sc%rep(i)%labels%type = islabels
                         call w%sc%rep(i)%labels%style%reset(w%sc%rep(i))
                      end if
                   end if
                elseif (w%sc%rep(i)%type == reptype_unitcell) then
                   if (changedisplay(4)) w%sc%rep(i)%shown = isuc
                end if
             end if
          end do
          chbuild = .true.
       end if
    end if

    ! scene menu
    str1="##viewscenebutton" // c_null_char
    if (iw_button("Scene",disabled=.not.associated(w%sc))) then
       call igOpenPopup_Str(c_loc(str1),ImGuiPopupFlags_None)
    end if
    if (associated(w%sc)) then
       if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_None)) then
          ! display shortcuts
          call iw_text("Display Shortcuts",highlight=.true.)
          if (iw_checkbox("Atoms##atomsshortcut",isatom)) then
             do i = 1, w%sc%nrep
                if (w%sc%rep(i)%isinit) then
                   if (w%sc%rep(i)%type == reptype_atoms) then
                      w%sc%rep(i)%atoms%display = isatom
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display atoms in all objects ("//&
             trim(get_bind_keyname(BIND_VIEW_TOGGLE_ATOMS)) // ").",ttshown)

          if (iw_checkbox("Bonds##bondsshortcut",isbond,sameline=.true.)) then
             do i = 1, w%sc%nrep
                if (w%sc%rep(i)%isinit) then
                   if (w%sc%rep(i)%type == reptype_atoms) then
                      w%sc%rep(i)%bonds%display = isbond
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display bonds in all objects ("//&
             trim(get_bind_keyname(BIND_VIEW_TOGGLE_BONDS)) // ").",ttshown)

          if (iw_checkbox("Labels##labelshortcut",islabelsl,sameline=.true.)) then
             do i = 1, w%sc%nrep
                if (w%sc%rep(i)%isinit) then
                   if (w%sc%rep(i)%type == reptype_atoms) then
                      w%sc%rep(i)%labels%display = islabelsl
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display labels in all objects ("//&
             trim(get_bind_keyname(BIND_VIEW_CYCLE_LABELS)) // ").",ttshown)

          if (.not.sys(w%view_selected)%c%ismolecule) then
             if (iw_checkbox("Unit Cell##ucshortcut",isuc,sameline=.true.)) then
                do i = 1, w%sc%nrep
                   if (w%sc%rep(i)%isinit) then
                      if (w%sc%rep(i)%type == reptype_unitcell) then
                         w%sc%rep(i)%shown = isuc
                      end if
                   end if
                end do
                chbuild = .true.
             end if
             call iw_tooltip("Toggle display unit cell in all objects ("//&
                trim(get_bind_keyname(BIND_VIEW_TOGGLE_CELL)) // ").",ttshown)
          end if

          if (iw_checkbox("Polyhedra##polyshortcut",ispoly,sameline=.true.)) then
             do i = 1, w%sc%nrep
                if (w%sc%rep(i)%isinit) then
                   if (w%sc%rep(i)%type == reptype_atoms) then
                      w%sc%rep(i)%poly%display = ispoly
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display polyhedra in all objects ("//&
             trim(get_bind_keyname(BIND_VIEW_TOGGLE_POLYHEDRA)) // ").",ttshown)

          ! periodicity (number of cells) selector
          if (.not.sys(w%view_selected)%c%ismolecule) then
             ! title
             call igAlignTextToFramePadding()
             call iw_text("Periodicity",highlight=.true.)
             if (iw_button("Reset##periodicity",sameline=.true.)) then
                w%sc%nc = 1
                chbuild = .true.
             end if
             call iw_tooltip("Reset the number of cells to one",ttshown)

             ! number of cells in each direction (shared digit width)
             nc = w%sc%nc
             ldum = iw_intstepper("aaxis",nc(1),label="a:",minval=1_c_int,notlive=.true.,&
                tooltip="Number of cells represented along the a crystallographic axis")
             ldum = iw_intstepper("baxis",nc(2),label="b:",minval=1_c_int,sameline=.true.,notlive=.true.,&
                tooltip="Number of cells represented along the b crystallographic axis")
             ldum = iw_intstepper("caxis",nc(3),label="c:",minval=1_c_int,sameline=.true.,notlive=.true.,&
                tooltip="Number of cells represented along the c crystallographic axis")
             if (any(nc /= w%sc%nc)) then
                w%sc%nc = nc
                chbuild = .true.
             end if
          end if

          ! camera position: align view axis
          call iw_text("Camera Position",highlight=.true.)
          if (.not.sys(w%view_selected)%c%ismolecule) then
             if (iw_button("a")) then
                call w%sc%align_view_axis(1)
                chrender = .true.
             end if
             call iw_tooltip("Align the camera along the crystallographic a axis ("//&
                trim(get_bind_keyname(BIND_VIEW_ALIGN_A_AXIS)) // ").",ttshown)
             if (iw_button("b",sameline=.true.)) then
                call w%sc%align_view_axis(2)
                chrender = .true.
             end if
             call iw_tooltip("Align the camera along the crystallographic b axis ("//&
                trim(get_bind_keyname(BIND_VIEW_ALIGN_B_AXIS)) // ").",ttshown)
             if (iw_button("c",sameline=.true.)) then
                call w%sc%align_view_axis(3)
                chrender = .true.
             end if
             call iw_tooltip("Align the camera along the crystallographic c axis ("//&
                trim(get_bind_keyname(BIND_VIEW_ALIGN_C_AXIS)) // ").",ttshown)
          end if
          if (iw_button("x",sameline=.not.sys(w%view_selected)%c%ismolecule)) then
             call w%sc%align_view_axis(-1)
             chrender = .true.
          end if
          call iw_tooltip("Align the camera along the Cartesian x axis ("//&
             trim(get_bind_keyname(BIND_VIEW_ALIGN_X_AXIS)) // ").",ttshown)
          if (iw_button("y",sameline=.true.)) then
             call w%sc%align_view_axis(-2)
             chrender = .true.
          end if
          call iw_tooltip("Align the camera along the Cartesian y axis ("//&
             trim(get_bind_keyname(BIND_VIEW_ALIGN_Y_AXIS)) // ").",ttshown)
          if (iw_button("z",sameline=.true.)) then
             call w%sc%align_view_axis(-3)
             chrender = .true.
          end if
          call iw_tooltip("Align the camera along the Cartesian z axis ("//&
             trim(get_bind_keyname(BIND_VIEW_ALIGN_Z_AXIS)) // ").",ttshown)
          ch = iw_dragfloat_realc("Reset Distance##resetdistance",x1=w%sc%camresetdist,speed=0.01_c_float,&
             min=0.1_c_float,max=8.0_c_float,decimal=2,sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
          if (ch) chrender = .true. ! constant-size labels and the gizmo depend on camresetdist
          call iw_tooltip("Ratio controlling distance from object when resetting camera",ttshown)

          ! projection mode
          if (iw_radiobutton("Orthographic##projortho",bool=w%sc%isortho,boolval=.true.)) then
             call w%sc%update_projection_matrix()
             chrender = .true.
          end if
          call iw_tooltip("Parallel projection, with no perspective distortion",ttshown)
          if (iw_radiobutton("Perspective##projpersp",bool=w%sc%isortho,boolval=.false.,sameline=.true.)) then
             call w%sc%update_projection_matrix()
             chrender = .true.
          end if
          call iw_tooltip("Perspective projection, with distant objects appearing smaller",ttshown)

          ! object resolution
          call iw_text("Object Resolution",highlight=.true.)
          ires = w%sc%atom_res - 1
          call iw_combo_simple("Atoms##atomresselect","1: Carnby"//c_null_char//"2: Rough"//c_null_char//&
             "3: Normal"//c_null_char//"4: Good"//c_null_char//"5: Amazing"//c_nulL_char,ires)
          call iw_tooltip("Set the resolution of the spheres representing the atoms",ttshown)
          if (ires + 1 /= w%sc%atom_res) then
             w%sc%atom_res = ires + 1
             chrender = .true.
          end if
          ires = w%sc%bond_res - 1
          call iw_combo_simple("Bonds##bondresselect","1: Rough" // c_null_char // "2: Normal" // c_null_char //&
             "3: Good" // c_null_char,ires,sameline=.true.)
          call iw_tooltip("Set the resolution of the cylinders representing the bonds",ttshown)
          if (ires + 1 /= w%sc%bond_res) then
             w%sc%bond_res = ires + 1
             chrender = .true.
          end if

          ! scene appearance (atom border color is set per representation in
          ! the edit-representations window)
          call iw_text("Appearance",highlight=.true.)

          ! background color
          chrender = chrender .or. iw_coloredit("Background",rgb=w%sc%bgcolor)
          call iw_tooltip("Change the scene background color",ttshown)

          ! apply to all scenes
          if (iw_button("Apply to All Systems",danger=.true.)) then
             do i = 1, nsys
                if (sysc(i)%status == sys_init .and. i /= w%view_selected) then
                   ! atoms, bonds, unit cell
                   do j = 1, sysc(i)%sc%nrep
                      if (sysc(i)%sc%rep(j)%isinit) then
                         if (sysc(i)%sc%rep(j)%type == reptype_atoms) then
                            sysc(i)%sc%rep(j)%atoms%display = isatom
                            sysc(i)%sc%rep(j)%bonds%display = isbond
                            sysc(i)%sc%rep(j)%poly%display = ispoly
                            sysc(i)%sc%rep(j)%labels%display = islabelsl
                            if (islabelsl) then
                               if (sys(i)%c%ismolecule.and.islabels == 8) then
                                  sysc(i)%sc%rep(j)%labels%type = 0
                               else
                                  sysc(i)%sc%rep(j)%labels%type = islabels
                               end if
                               call sysc(i)%sc%rep(j)%labels%style%reset(sysc(i)%sc%rep(j))
                            end if
                         elseif (sysc(i)%sc%rep(j)%type == reptype_unitcell.and.&
                            .not.sys(w%view_selected)%c%ismolecule) then
                            sysc(i)%sc%rep(j)%shown = isuc
                         end if
                      end if
                   end do
                   ! rest
                   if (.not.sys(w%view_selected)%c%ismolecule.and..not.sys(i)%c%ismolecule) &
                      sysc(i)%sc%nc = w%sc%nc
                   sysc(i)%sc%atom_res = w%sc%atom_res
                   sysc(i)%sc%bond_res = w%sc%bond_res
                   sysc(i)%sc%bgcolor = w%sc%bgcolor
                   sysc(i)%sc%camresetdist = w%sc%camresetdist
                   sysc(i)%sc%isortho = w%sc%isortho
                   if (sysc(i)%sc%iscaminit) call sysc(i)%sc%update_projection_matrix()
                end if
                call sysc(i)%sc%build_lists()
             end do
          end if
          call iw_tooltip("Apply these settings to all systems",ttshown)
          if (iw_button("Reset",sameline=.true.,danger=.true.)) then
             call w%sc%init(w%view_selected)
             chbuild = .true.
          end if
          call iw_tooltip("Reset to the default settings",ttshown)

          call igEndPopup()
       end if
    end if
    call iw_tooltip("Change the view options",ttshown)

    ! gear menu
    str1="##viewgear" // c_null_char
    if (iw_button("Objects",sameline=.true.,disabled=.not.associated(w%sc))) then
       call igOpenPopup_Str(c_loc(str1),ImGuiPopupFlags_None)
    end if
    if (associated(w%sc)) then
       if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_None)) then
          call igAlignTextToFramePadding()
          ! objects table
          call iw_text("List of Objects",highlight=.true.)

          ! add button
          ldum = iw_button("Add",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
          if (ok) then
             if (iw_menuitem("Ball and Stick",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_ballandstick)
                chbuild = .true.
             end if
             call iw_tooltip("Draw atoms as balls and bonds as sticks, hide the labels",ttshown)

             if (iw_menuitem("Bonds",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_sticks)
                chbuild = .true.
             end if
             call iw_tooltip("Draw bonds as sticks, hide atoms and labels",ttshown)

             if (iw_menuitem("Licorice",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_licorice)
                chbuild = .true.
             end if
             call iw_tooltip("Draw atoms and bonds with the same radius, hide labels",ttshown)

             call igSeparator()
             if (iw_menuitem("Van der Waals Contacts",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_vdwcontacts)
                chbuild = .true.
             end if
             call iw_tooltip("Display contacts between nonbonded atoms closer than the sum &
                &of their van der Waals radii",ttshown)

             if (iw_menuitem("Hydrogen Bonds",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_hbonds)
                chbuild = .true.
             end if
             call iw_tooltip("Display contacts between hydrogen bonded atoms",ttshown)

             call igSeparator()
             if (iw_menuitem("Critical Points",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_criticalpoints)
                chbuild = .true.
             end if
             call iw_tooltip("Draw dummy atoms representing critical points (Xn, Xb,... atoms)",ttshown)

             if (iw_menuitem("Gradient Paths",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_gradientpaths)
                chbuild = .true.
             end if
             call iw_tooltip("Draw dummy atoms representing gradient paths (Xz atoms)",ttshown)

             call igSeparator()
             if (iw_menuitem("Coordination Polyhedra",shortcut_text="Atoms")) then
                call w%sc%add_representation(reptype_atoms,repflavor_atoms_polyhedra)
                chbuild = .true.
             end if
             call iw_tooltip("Draw coordination polyhedra around the center atoms",ttshown)

             if (.not.sys(w%view_selected)%c%ismolecule) then
                call igSeparator()
                if (iw_menuitem("Unit Cell",shortcut_text="Cell")) then
                   call w%sc%add_representation(reptype_unitcell,repflavor_unitcell_basic)
                   chbuild = .true.
                end if
                call iw_tooltip("Display the unit cell",ttshown)
             end if

             call igSeparator()
             if (iw_menuitem("Cartesian/Crystallographic Axes",shortcut_text="Axes")) then
                call w%sc%add_representation(reptype_axes,repflavor_axes)
                chbuild = .true.
             end if
             call iw_tooltip("Display the cartesian (lab-frame) x/y/z axes",ttshown)

             ! symmetry available for crystals (always) or molecules with a point group
             call igSeparator()
             enabled = .true.
             if (sys(w%view_selected)%c%ismolecule) enabled = sys(w%view_selected)%c%pg%avail
             if (iw_menuitem("Symmetry Elements",shortcut_text="Symmetry",enabled=enabled)) then
                call w%sc%add_representation(reptype_symelem,repflavor_symelem)
                chbuild = .true.
             end if
             call iw_tooltip("Display symmetry elements",ttshown)

             call igEndPopup()
          end if
          call iw_tooltip("Add a new object to the view",ttshown)

          ! set table style
          sz%x = 3._c_float
          sz%y = 1._c_float
          call igPushStyleVar_Vec2(ImGuiStyleVar_FramePadding,sz)
          sz%x = 8._c_float
          sz%y = 5._c_float
          call igPushStyleVar_Vec2(ImGuiStyleVar_ItemSpacing,sz)
          sz%x = 2._c_float
          sz%y = 2._c_float
          call igPushStyleVar_Vec2(ImGuiStyleVar_CellPadding,sz)

          ! rest of the table
          str2 = "Objects##0,0" // c_null_char
          flags = ImGuiTableFlags_NoSavedSettings
          flags = ior(flags,ImGuiTableFlags_RowBg)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          sz0%x = 0
          nrep = count(w%sc%rep(1:w%sc%nrep)%isinit)
          nrep = min(nrep,10)
          sz0%y = iw_calcheight(nrep,0,.true.)
          if (igBeginTable(c_loc(str2),5,flags,sz0,0._c_float)) then
             str3 = "##1closebutton" // c_null_char
             flags = ImGuiTableColumnFlags_None
             width = max(4._c_float, fontsize%y + 2._c_float)
             call igTableSetupColumn(c_loc(str3),flags,width,ic_closebutton)

             str3 = "##1viewbutton" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str3),flags,0.0_c_float,ic_viewbutton)

             str3 = "Name##1name" // c_null_char
             flags = ImGuiTableColumnFlags_WidthStretch
             call igTableSetupColumn(c_loc(str3),flags,0.0_c_float,ic_name)

             str3 = "Type##1type" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str3),flags,0.0_c_float,ic_type)

             str3 = "##1editbutton" // c_null_char
             flags = ImGuiTableColumnFlags_None
             width = iw_calcwidth(4,1)
             call igTableSetupColumn(c_loc(str3),flags,width,ic_editbutton)

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             if (w%sc%representation_menu(w%id)) chbuild = .true.

             call igEndTable()
          end if
          call igPopStyleVar(3_c_int)

          call igEndPopup()
       end if
    end if
    call iw_tooltip("Add, remove, and modify objects",ttshown)

    ! update the draw lists and render
    if (associated(w%sc)) then
       if (w%sc%timelastbuild < sysc(w%view_selected)%timelastchange_buildlists) then
          w%sc%forcebuildlists = .true.
          ! during interactive dynamics the geometry changes every frame but the
          ! atom identities are stable, so keep the pick index (else grabbing an
          ! atom would be impossible)
          if (.not.sysc(w%view_selected)%md_run) w%mousepos_idx = 0
       end if
       if (chbuild) w%sc%forcebuildlists = .true.
       if (chrender .or. w%sc%forcebuildlists .or. w%sc%timelastrender < sysc(w%view_selected)%timelastchange_render) &
          w%forcerender = .true.

       ! continuous render if animation is active
       if (w%sc%ifreq_selected > 0.and.w%sc%iqpt_selected > 0.and.sys(w%view_selected)%c%vib%hasvibs.and.&
          w%sc%animation > 0) &
          w%forcerender = .true.
    end if

    ! export image
    ldum = iw_button("Tools",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft,&
       disabled=.not.associated(w%sc))
    call iw_tooltip("Show various tools operating on the view of this system",ttshown)
    if (ok) then
       enabled = associated(w%sc)
       if (iw_menuitem("Export to Image...",enabled=enabled)) &
          iaux = stack_create_window(wintype_exportimage,.true.,idparent=w%id,orraise=-1)
       call iw_tooltip("Export the current view to an image file (png)",ttshown)

       ! separator
       call igSeparator()

       if (iw_menuitem("View/Edit Geometry...",enabled=enabled)) &
          iaux = stack_create_window(wintype_geometry,.true.,isys=w%view_selected,orraise=-1)
       call iw_tooltip("View and edit the atomic positions, bonds, etc.",ttshown)

       if (iw_menuitem("Recalculate bonds",BIND_RECALC_BONDS,enabled=enabled)) &
          call sysc(w%view_selected)%rebond()
       call iw_tooltip("Recompute the bonds/connectivity for this system",ttshown)

       if (iw_menuitem("Vibrations...",enabled=enabled)) &
          iaux = stack_create_window(wintype_vibrations,.true.,idparent=w%id,orraise=-1)
       call iw_tooltip("Display an animation showing the atomic vibrations for this system",ttshown)

       if (iw_menuitem("Dynamics...",enabled=enabled)) &
          iaux = stack_create_window(wintype_dynamics,.true.,idparent=w%id,orraise=-1)
       call iw_tooltip("Run an interactive molecular-dynamics simulation: animate the system at a &
          &given temperature and drag atoms with the mouse",ttshown)
       call igEndPopup()
    end if

    ! camera lock
    if (w%ismain) then
       ldum = iw_button("Cam-Lock",disabled=.not.associated(w%sc),sameline=.true.,&
          popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
       if (ok) then
          if (iw_menuitem("Lock All",selected=(lockbehavior==2))) then
             lockbehavior = 2
             do k = 1, nsys
                sysc(k)%sc%lockedcam = -1
             end do
          end if
          call iw_tooltip("Lock the camera position for all loaded systems",ttshown)

          if (iw_menuitem("Lock SCF Iterations Only",selected=(lockbehavior==1))) then
             lockbehavior = 1
             do k = 1, nsys
                if (sysc(k)%collapse < 0) then
                   sysc(k)%sc%lockedcam = k
                elseif (sysc(k)%collapse > 0) then
                   sysc(k)%sc%lockedcam = sysc(k)%collapse
                else
                   sysc(k)%sc%lockedcam = 0
                end if
             end do
          end if
          call iw_tooltip("Lock the camera position only for SCF iterations of the same system",ttshown)

          if (iw_menuitem("Unlock All",selected=(lockbehavior==0))) then
             lockbehavior = 0
             do k = 1, nsys
                sysc(k)%sc%lockedcam = 0
             end do
          end if
          call iw_tooltip("Do not lock the camera position for any system",ttshown)

          call igEndPopup()
       end if
       call iw_tooltip("Lock the camera position and orientation for multiple systems",ttshown)
    end if

    ! the button for new alternate view
    if (iw_button("+",disabled=.not.associated(w%sc),sameline=.true.)) then
       idum = stack_create_window(wintype_view,.true.,purpose=wpurp_view_alternate)
       win(idum)%sc = w%sc
       ! the value copy aliased the source scene's GL handles; detach so the new
       ! view lazily builds its own instance buffers
       call win(idum)%sc%gl%detach()
       win(idum)%view_selected = w%view_selected
       call win(idum)%sc%reset_animation()
    end if
    call iw_tooltip("Create a new view for the current scene",ttshown)

    ! the selected system combo
    call igSameLine(0._c_float,-1._c_float)
    str2 = "" // c_null_char
    if (goodsys) then
       str2 = string(w%view_selected) // ": " // trim(sysc(w%view_selected)%seed%name) // c_null_char
    end if
    str1 = "##systemcombo" // c_null_char
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nsys
          if (sysc(i)%status == sys_init) then
             is_selected = (w%view_selected == i)
             str2 = string(i) // ": " // trim(sysc(i)%seed%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) &
                call w%select_view(i)
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Choose the system displayed",ttshown)

    ! get the remaining size for the texture
    call igGetContentRegionAvail(szavail)
    szavail%y = szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y

    ! resize the render texture to the window (snapped to a bucket and capped),
    ! so the still image stays sharp without pixelation on large windows. Only
    ! reallocate when the bucketed side changes, to avoid churn during resize.
    newside = max(ceiling(max(szavail%x,szavail%y)),1)
    newside = ((newside + texture_side_bucket - 1) / texture_side_bucket) * texture_side_bucket
    newside = min(max(newside, min_texture_side), max_texture_side)
    if (newside /= w%FBOside) then
       call w%delete_texture_view()
       call w%create_texture_view(newside)
       w%forcerender = .true.
    end if

    ! Adaptive resolution: while the camera is being manipulated, render at a
    ! reduced (interactive) resolution to stay responsive on large windows. When
    ! the user stops, snap back to a full-resolution render. Whenever the desired
    ! render scale (interacting) differs from what the texture currently holds
    ! (lowresrender), force a render this frame so the displayed sub-region (the
    ! UVs below, scaled by rscale) always matches the texture content. Otherwise
    ! a press that switches to low-res before any motion would sample the
    ! sub-region of a full-region render (scene out of position).
    interacting = associated(w%sc) .and. (w%ilock /= ilock_no)
    if (interacting .neqv. w%lowresrender) w%forcerender = .true.
    if (interacting) then
       rscale = real(min(w%FBOside,interactive_texture_side),c_float) / real(w%FBOside,c_float)
    else
       rscale = 1._c_float
    end if

    ! draw the texture, largest region with the same shape as the available region
    ! that fits into the texture square
    scal = real(w%FBOside,c_float) / max(max(szavail%x,szavail%y),1._c_float)
    sz0%x = 0.5 * (real(w%FBOside,c_float) - szavail%x * scal) / real(w%FBOside,c_float)
    sz0%y = 0.5 * (real(w%FBOside,c_float) - szavail%y * scal) / real(w%FBOside,c_float)
    sz1%x = 1._c_float - sz0%x
    sz1%y = 1._c_float - sz0%y

    ! record the visible (cropped) region so the scene can place
    ! window-anchored objects (e.g. the axes gizmo) relative to it. This uses the
    ! unscaled crop (relative to the full square render). If the visible region
    ! changed and the scene has window-anchored objects, force a re-render so
    ! they track the new window geometry.
    if (associated(w%sc)) then
       if (w%sc%hasanchoredobj .and. &
          (abs(sz0%x - w%sc%viewuv0(1)) > 1e-6_c_float .or. abs(sz0%y - w%sc%viewuv0(2)) > 1e-6_c_float)) &
          w%forcerender = .true.
       w%sc%viewuv0 = (/sz0%x,sz0%y/)
    end if

    ! during interactive low-res rendering the scene is rasterized into the
    ! lower-left rscale-fraction of the texture, so sample only that sub-region.
    sz0%x = rscale * sz0%x
    sz0%y = rscale * sz0%y
    sz1%x = rscale * sz1%x
    sz1%y = rscale * sz1%y

    ! render the image to the texture, if requested
    if (w%forcerender) then
       ! viewport side: full texture, or the interactive sub-square (lower-left)
       vside = max(nint(rscale * w%FBOside),1)

       ! render to the draw framebuffer
       call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
       call glViewport(0_c_int,0_c_int,vside,vside)
       if (associated(w%sc)) then
          call glClearColor(w%sc%bgcolor(1),w%sc%bgcolor(2),&
             w%sc%bgcolor(3),1._c_float)
       else
          call glClearColor(ColorClearTransparent(1),ColorClearTransparent(2),&
             ColorClearTransparent(3),ColorClearTransparent(4))
       end if
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
       if (associated(w%sc)) call w%sc%render()
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! render to the pick frame buffer. Skip this (heavy RGBA32F) pass while
       ! interacting: picking is only queried when the view is idle, and the
       ! mouse-to-texture mapping assumes the full-resolution pick buffer.
       if (.not.interacting) then
          call glBindFramebuffer(GL_FRAMEBUFFER, w%FBOpick)
          call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
          call glClearColor(ColorClearTransparent(1),ColorClearTransparent(2),&
             ColorClearTransparent(3),ColorClearTransparent(4))
          call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
          if (associated(w%sc)) call w%sc%renderpick()
          call glBindFramebuffer(GL_FRAMEBUFFER, 0)
       end if

       w%lowresrender = interacting
       w%forcerender = .false.
    end if

    ! border and tint for the image, draw the image, update the rectangle
    tintcol = ColorWhite
    bgcol = ColorBlack
    call igPushStyleColor_Vec4(ImGuiCol_Button,bgcol)
    call igPushStyleColor_Vec4(ImGuiCol_ButtonActive,bgcol)
    call igPushStyleColor_Vec4(ImGuiCol_ButtonHovered,bgcol)
    ! the render texture is a bottom-left-origin OpenGL FBO; flip the vertical
    ! UVs so the scene is presented right-side up. Without this the 3D geometry
    ! (and the Cartesian-axes gizmo) renders vertically mirrored / left-handed.
    tmpuv = sz0%y
    sz0%y = sz1%y
    sz1%y = tmpuv
    str1 = "##imagebutton" // c_null_char
    ldum = igImageButtonEx(igGetID_Str(c_loc(str1)),w%FBOtex, szavail, sz0, sz1, szero, bgcol, tintcol)
    call igPopStyleColor(3)

    ! get view geometry and mouse position
    call igGetItemRectMin(w%v_rmin)
    call igGetItemRectMax(w%v_rmax)
    call igGetMousePos(pos)

    ! set view modes based on user's key presses
    if (goodsys) call w%viewmode_set_mode()

    ! get hover and needpick
    hover = goodsys .and. w%ilock == ilock_no
    if (hover) hover = igIsItemHovered(ImGuiHoveredFlags_None)
    needpick = hover .and. (abs(w%mousepos_lastpick%x-pos%x) > 1e-4.or.abs(w%mousepos_lastpick%y-pos%y) > 1e-4) .and.&
       (w%timelast_view_getpixel + pick_interval < time)
    needpick = needpick .or. w%viewmode_activate_picking(hover)
    ! during interactive dynamics, re-pick each throttle interval even if the
    ! cursor is still, so the atom under a stationary cursor stays selectable as
    ! it moves (needed for grab-and-drag)
    if (hover .and. w%view_selected >= 1 .and. w%view_selected <= nsys) needpick = needpick .or. &
       (sysc(w%view_selected)%md_run .and. (w%timelast_view_getpixel + pick_interval < time))

    ! get the ID of the atom under mouse
    if (needpick) then
       w%mousepos_idx = 0
       w%mousepos_lastpick = pos
       call w%mousepos_to_texpos(pos)
       call w%getpixel(pos,rgba=rgba)

       ! transform to atom cell ID and lattice vector
       w%mousepos_idx(1:4) = transfer(rgba,w%mousepos_idx(1:4))
       w%mousepos_idx(5) = w%mousepos_idx(1)
       if (associated(w%sc) .and. w%mousepos_idx(1) > 0 .and. w%mousepos_idx(1) <= w%sc%obj%nsph) then
          w%mousepos_idx(1:4) = w%sc%obj%sph(w%mousepos_idx(1))%idx
       else
          w%mousepos_idx = 0
       end if
       w%timelast_view_getpixel = time
    elseif (.not.hover) then
       w%mousepos_idx = 0
    end if

    ! the viewmode display on the bar
    call w%viewmode_bar_display()

    ! atom hover message
    if (hover .and. w%mousepos_idx(1) > 0) then
       call igSameLine(0._c_float,-1._c_float)
       icel = w%mousepos_idx(1)
       is = sys(w%view_selected)%c%atcel(icel)%is
       ineq = sys(w%view_selected)%c%atcel(icel)%idx
       ismol = sys(w%view_selected)%c%ismolecule

       ! lead with the occupant list for a mixed site, otherwise the atom name
       ! (mix_string returns empty for a single-occupant site)
       msg = sys(w%view_selected)%c%mix_string(ineq)
       if (len_trim(msg) == 0) &
          msg = trim(sys(w%view_selected)%c%at(ineq)%name)
       if (.not.ismol) then
          x0 = sys(w%view_selected)%c%atcel(icel)%x

          msg = trim(msg) // " [cellid=" // string(icel) // "+(" // string(w%mousepos_idx(2)) // "," // string(w%mousepos_idx(3)) //&
             "," // string(w%mousepos_idx(4)) // "),nneqid=" // string(ineq) // ",wyckoff=" // &
             string(sys(w%view_selected)%c%at(ineq)%mult) // string(sys(w%view_selected)%c%at(ineq)%wyc)
          if (sys(w%view_selected)%c%nmol > 1) &
             msg = msg // ",molid=" // string(sys(w%view_selected)%c%idatcelmol(1,icel))
          msg = msg // "] " //&
             string(x0(1)+w%mousepos_idx(2),'f',decimal=4) //" "// string(x0(2)+w%mousepos_idx(3),'f',decimal=4) //" "//&
             string(x0(3)+w%mousepos_idx(4),'f',decimal=4) // " (frac)"
       else
          x0 = (sys(w%view_selected)%c%atcel(icel)%r+sys(w%view_selected)%c%molx0) * dunit0(iunit_ang)

          msg = trim(msg) // " [id=" // string(icel)
          if (sys(w%view_selected)%c%nmol > 1) &
             msg = msg // ",molid=" // string(sys(w%view_selected)%c%idatcelmol(1,icel))
          msg = msg // "] " //&
             string(x0(1),'f',decimal=4) //" "// string(x0(2),'f',decimal=4) //" "//&
             string(x0(3),'f',decimal=4) // " (Å)"
       end if
       call iw_text(msg)
    end if

    ! tooltip for distance measurement
    if (hover) &
       call w%draw_selection_tooltip(w%mousepos_idx)

    ! Process mouse events
    call w%viewmode_process_events(hover)

    ! process keybindings
    !! increase and decrease the number of cells in main view
    if (associated(w%sc)) then
       if (w%ismain) then
          if (.not.sys(w%view_selected)%c%ismolecule) then
             if (is_bind_event(BIND_VIEW_INC_NCELL)) then
                do i = 1, 3
                   if (sys(w%view_selected)%c%vaclength(i) < iperiod_vacthr) &
                      w%sc%nc(i) = w%sc%nc(i) + 1
                end do
                w%sc%forcebuildlists = .true.
             elseif (is_bind_event(BIND_VIEW_DEC_NCELL)) then
                do i = 1, 3
                   if (sys(w%view_selected)%c%vaclength(i) < iperiod_vacthr) &
                      w%sc%nc(i) = w%sc%nc(i) - 1
                end do
                w%sc%nc = max(w%sc%nc,1)
                w%sc%forcebuildlists = .true.
             end if
             if (is_bind_event(BIND_VIEW_ALIGN_A_AXIS)) then
                call w%sc%align_view_axis(1)
                w%forcerender = .true.
             end if
             if (is_bind_event(BIND_VIEW_ALIGN_B_AXIS)) then
                call w%sc%align_view_axis(2)
                w%forcerender = .true.
             end if
             if (is_bind_event(BIND_VIEW_ALIGN_C_AXIS)) then
                call w%sc%align_view_axis(3)
                w%forcerender = .true.
             end if
          end if
          if (is_bind_event(BIND_VIEW_ALIGN_X_AXIS)) then
             call w%sc%align_view_axis(-1)
             w%forcerender = .true.
          end if
          if (is_bind_event(BIND_VIEW_ALIGN_Y_AXIS)) then
             call w%sc%align_view_axis(-2)
             w%forcerender = .true.
          end if
          if (is_bind_event(BIND_VIEW_ALIGN_Z_AXIS)) then
             call w%sc%align_view_axis(-3)
             w%forcerender = .true.
          end if
       end if
    end if

    ! right-align the rest of the contents
    call igSameLine(0._c_float,0._c_float)
    call iw_setposx_fromend(5,1)

    ! keyboard actions on the current atom selection, when the view is focused:
    selcleared = .false.
    if (w%focused() .and. ok_system(w%view_selected,sys_init)) then
       is = w%view_selected
       if (allocated(sysc(is)%highlight_rgba)) then
          ok = any(sysc(is)%highlight_rgba >= 0._c_float)
       else
          ok = .false.
       end if
       if (ok .and. is_bind_event(BIND_EDITSELECT_REMOVE)) then
          ! delete the selected atoms
          call sysc(is)%edit_highlighted_atoms(remove=.true.,errmsg=msg)
          sysc(is)%sc%nextbuildlists_fixcam = .true.
          w%forcerender = .true.
       elseif (ok .and. is_bind_event(BIND_EDITSELECT_DESELECT)) then
          ! clear the selection
          call sysc(is)%highlight_clear(.false.)
          w%forcerender = .true.
          selcleared = .true.
       elseif (is_bind_event(BIND_EDITSELECT_SELECT_ALL)) then
          ! select all atoms
          call sysc(is)%highlight_all()
          w%forcerender = .true.
       end if
    end if

    if (.not.w%ismain) then
       ! the close button
       if (iw_button("Close",danger=.true.)) w%isopen = .false.

       ! exit if focused and received the close keybinding (unless Escape was
       ! just used to clear the selection)
       if (.not.selcleared) then
          if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) .or.&
             is_bind_event(BIND_CLOSE_ALL_DIALOGS)) then
             w%isopen = .false.
          end if
       end if
    end if

  end subroutine draw_view

  !> Create the texture for the view window, with atex x atex pixels.
  module subroutine create_texture_view(w,atex)
    use interfaces_opengl3
    use gui_main, only: ColorClearTransparent
    use tools_io, only: ferror, faterr
    class(window), intent(inout), target :: w
    integer, intent(in) :: atex

    ! FBO and buffers
    call glGenTextures(1, c_loc(w%FBOtex))
    call glGenTextures(1, c_loc(w%FBOrgba))
    call glGenRenderbuffers(1, c_loc(w%FBOdepth))
    call glGenRenderbuffers(1, c_loc(w%FBOdepthp))
    call glGenFramebuffers(1, c_loc(w%FBO))
    call glGenFramebuffers(1, c_loc(w%FBOpick))

    ! textures
    call glBindTexture(GL_TEXTURE_2D, w%FBOtex)
    call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, atex, atex, 0, GL_RGBA, GL_UNSIGNED_BYTE, c_null_ptr)
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    call glBindTexture(GL_TEXTURE_2D, 0)

    call glBindTexture(GL_TEXTURE_2D, w%FBOrgba)
    call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, atex, atex, 0, GL_RGBA, GL_FLOAT, c_null_ptr)
    call glBindTexture(GL_TEXTURE_2D, 0)

    ! render buffers
    call glBindRenderbuffer(GL_RENDERBUFFER, w%FBOdepth)
    call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, atex, atex)
    call glBindRenderbuffer(GL_RENDERBUFFER, 0)
    call glBindRenderbuffer(GL_RENDERBUFFER, w%FBOdepthp)
    call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, atex, atex)
    call glBindRenderbuffer(GL_RENDERBUFFER, 0)

    ! frame buffers
    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
    call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, w%FBOtex, 0)
    call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, w%FBOdepth)
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
       call ferror('window_init','framebuffer (draw) is not complete',faterr)
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBOpick)
    call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, w%FBOrgba, 0)
    call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, w%FBOdepthp)
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
       call ferror('window_init','framebuffer (pick) is not complete',faterr)
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

    ! write the texture side
    w%FBOside = atex

    ! initial clear
    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
    call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
    call glClearColor(ColorClearTransparent(1),ColorClearTransparent(2),&
       ColorClearTransparent(3),ColorClearTransparent(4))
    call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBOpick)
    call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
    call glClearColor(ColorClearTransparent(1),ColorClearTransparent(2),&
       ColorClearTransparent(3),ColorClearTransparent(4))
    call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

  end subroutine create_texture_view

  !> Delete the texture for the view window
  module subroutine delete_texture_view(w)
    use interfaces_opengl3
    class(window), intent(inout), target :: w

    call glDeleteTextures(1, c_loc(w%FBOtex))
    call glDeleteTextures(1, c_loc(w%FBOrgba))
    call glDeleteRenderbuffers(1, c_loc(w%FBOdepth))
    call glDeleteRenderbuffers(1, c_loc(w%FBOdepthp))
    call glDeleteFramebuffers(1, c_loc(w%FBO))
    call glDeleteFramebuffers(1, c_loc(w%FBOpick))

  end subroutine delete_texture_view

  !> Select system isys in view window.
  module subroutine select_view(w,isys)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: nsys, sysc, sys_init
    class(window), intent(inout), target :: w
    integer, intent(in) :: isys

    if (isys < 1 .or. isys > nsys) return
    if (sysc(isys)%status /= sys_init) w%forcerender = .true. ! for removing the last system in tree
    if (w%view_selected == isys) return

    ! select and render the new scene
    w%view_selected = isys
    if (w%ismain) then
       w%sc => sysc(w%view_selected)%sc
    else
       if (.not.associated(w%sc)) allocate(w%sc)
       call w%sc%end()
       call w%sc%init(w%view_selected)
    end if
    w%forcerender = .true.

    ! reset the mouse variables
    w%ilock = ilock_no
    w%mousepos_lastpick%x = 0._c_float
    w%mousepos_lastpick%y = 0._c_float
    w%mousepos_idx = 0

    ! reset the viewmodes
    w%viewmode = vm_navigate
    w%viewmode_transient = .false.

    ! set the time
    w%timelast_view_assign = glfwGetTime()

  end subroutine select_view

  !> Process the user keybindings that set the viewmode
  module subroutine viewmode_set_mode(w)
    use keybindings, only: is_bind_event, BIND_VIEWMODE_SELECT, BIND_VIEWMODE_MOVEMOL,&
       BIND_VIEWMODE_MOVEATOM
    use systems, only: nsys, sysc
    class(window), intent(inout), target :: w

    logical :: ok
    integer :: id

    ! while an interactive dynamics run is active on the viewed
    ! system, force the dedicated MD interaction mode and lock it: the
    ! user cannot switch modes until the run stops.
    if (w%view_selected >= 1 .and. w%view_selected <= nsys) then
       if (sysc(w%view_selected)%md_run) then
          w%viewmode = vm_mddrag
          w%viewmode_transient = .false.
          return
       end if
    end if
    if (w%viewmode == vm_mddrag) then
       ! the run stopped: leave the forced mode
       w%viewmode = vm_navigate
       w%viewmode_transient = .false.
    end if

    ! if the viewmode is forced by another window, check that window is still valid
    if (w%viewmode < 0) then
       id = w%vmdata%owner
       ok = (id > 0 .and. id <= nwin)
       if (ok) ok = win(id)%isinit
       if (.not.ok) then
          w%viewmode = vm_navigate
          w%viewmode_transient = .false.
          w%vmdata%owner = 0
       end if
    end if

    ! the window_forced view mode cannot be overridden by keybindings
    if (w%viewmode < 0) return

    ! select mode (shift)
    if (is_bind_event(BIND_VIEWMODE_SELECT,held=.true.)) then
       w%viewmode = vm_select
       w%viewmode_transient = .true.
    end if

    ! move-molecules mode (ctrl)
    if (is_bind_event(BIND_VIEWMODE_MOVEMOL,held=.true.)) then
       w%viewmode = vm_movemol
       w%viewmode_transient = .true.
    end if

    ! move-atoms mode (alt)
    if (is_bind_event(BIND_VIEWMODE_MOVEATOM,held=.true.)) then
       w%viewmode = vm_moveatom
       w%viewmode_transient = .true.
    end if

  end subroutine viewmode_set_mode

  !> Enter the window-forced transient view mode: the message is shown
  !> in the bar below the view and the mode only exits when the mouse
  !> clicks the view (any button) or a key is pressed (any key). If an atom
  !> is under the mouse when it is clicked, its ID is stored in
  !> w%viewmode_forced_idx (mousepos_idx layout) for use by the
  !> routine that set the mode; otherwise the result is zero. idcaller
  !> is the window ID (win(:)) of the caller, used by the caller to verify
  !> it owns the pick result.
  module subroutine viewmode_set_forced(w,mode,message,idcaller)
    class(window), intent(inout), target :: w
    integer, intent(in) :: mode
    character(len=*), intent(in) :: message
    integer, intent(in) :: idcaller

    w%viewmode = mode
    w%vmdata%owner = idcaller
    w%vmdata%msg = trim(message)
    w%vmdata%idx = 0

  end subroutine viewmode_set_forced

  !> Returns the tooltip message for the current viewmode
  module subroutine viewmode_bar_display(w)
    use gui_main, only: tooltip_delay, g
    use keybindings, only: get_bind_keyname, bindnames,&
       BIND_NUM, group_viewmode_navigation, group_viewmode_select,&
       group_viewmode_movemol, group_viewmode_moveatom, groupbind, BIND_NAV_ZOOM
    use utils, only: iw_combo_simple, iw_tooltip, igIsItemHovered_delayed, iw_text
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: ll, i, n, viewmode_before, iforced, nline
    character(len=:), allocatable :: viewmode_items
    integer, allocatable :: tips(:)
    character(len=64), allocatable :: keyline(:), lblline(:)
    logical :: ok

    logical, save :: ttshown = .false. ! tooltip flag

    viewmode_items = ""
    do i = 0, vm_NUM
       viewmode_items = trim(viewmode_items) // trim(vmnames(i)) // c_null_char
    end do

    if (w%viewmode < 0) then
       ! window-forced modes: show the combo; picking a normal mode
       ! leaves a window-forced pick, while vm_mddrag is re-forced
       ! each frame.
       iforced = vm_NUM+1
       viewmode_items = viewmode_items // trim(vmnames(w%viewmode)) // c_null_char
       call iw_combo_simple("##viewmode",viewmode_items,iforced)

       ! delayed tooltip describing the forced mode's mouse bindings
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
          call igBeginTooltip()
          if (w%viewmode == vm_mddrag) then
             call iw_text("Left drag",highlight=.true.)
             call iw_text("grab an atom / rotate the view",sameline=.true.)
             call iw_text("Right drag",highlight=.true.)
             call iw_text("move a molecule / translate the view",sameline=.true.)
             call iw_text("Middle drag",highlight=.true.)
             call iw_text("rotate a molecule",sameline=.true.)
          else
             call iw_text("Click",highlight=.true.)
             call iw_text("pick the atom under the cursor",sameline=.true.)
             call iw_text("Any key",highlight=.true.)
             call iw_text("cancel the pick",sameline=.true.)
          end if
          call igEndTooltip()
       end if

       if (iforced /= vm_NUM + 1) then
          w%viewmode = iforced
          w%viewmode_transient = .false.
       elseif (allocated(w%vmdata%msg)) then
          call iw_text(w%vmdata%msg,highlight=.true.,sameline=.true.)
       end if
    else
       ! the usual combo box; if the user selects the mode explicitly, it is not transient
       viewmode_before = w%viewmode
       call iw_combo_simple("##viewmode",viewmode_items,w%viewmode)
       if (w%viewmode /= viewmode_before) w%viewmode_transient = .false.

       ! get the tips for this view mode
       ll = 1
       n = 0
       allocate(tips(BIND_NUM))
       do i = 1, BIND_NUM
          ok = .false.
          if (w%viewmode == vm_navigate) then
             ok = (groupbind(i) == group_viewmode_navigation)
          elseif (w%viewmode == vm_select) then
             ok = (groupbind(i) == group_viewmode_select)
          elseif (w%viewmode == vm_movemol) then
             ok = (groupbind(i) == group_viewmode_movemol)
          elseif (w%viewmode == vm_moveatom) then
             ok = (groupbind(i) == group_viewmode_moveatom)
          end if
          if (ok) then
             n = n + 1
             tips(n) = i
             ll = max(ll,len_trim(get_bind_keyname(i)))
          end if
       end do

       ! delayed tooltip with info about the key/mouse bindings for this view mode
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) then
          if (igIsMouseHoveringRect(g%LastItemData%NavRect%min,g%LastItemData%NavRect%max,.false._c_bool)) then
             ! build the tooltip lines: one per bind in this mode's group, plus
             ! the mouse-scroll (cell-volume) line for the move modes
             allocate(keyline(n+1),lblline(n+1))
             do i = 1, n
                keyline(i) = trim(get_bind_keyname(tips(i)))
                lblline(i) = trim(bindnames(tips(i)))
             end do
             nline = n
             if (w%viewmode == vm_movemol .or. w%viewmode == vm_moveatom) then
                nline = nline + 1
                keyline(nline) = trim(get_bind_keyname(BIND_NAV_ZOOM))
                lblline(nline) = "Change cell volume (crystals)"
             end if
             ! align the key column
             ll = 1
             do i = 1, nline
                ll = max(ll,len_trim(keyline(i)))
             end do

             call igBeginTooltip()
             do i = 1, nline
                call iw_text(string(trim(keyline(i)),length=ll+1),highlight=.true.)
                call iw_text(trim(lblline(i)),sameline=.true.)
             end do
             call igEndTooltip()
             deallocate(keyline,lblline)
          end if
       end if
    end if

  end subroutine viewmode_bar_display

  !> Returns whether a pixel should be read from the picking baffer
  !> according to the current view mode and window state.  hover =
  !> whether the view is active and being hovered.
  module function viewmode_activate_picking(w,hover)
    use keybindings, only: is_bind_event, BIND_NAV_MEASURE, BIND_SELECT_MOLECULES_AND_DESELECT
    class(window), intent(inout), target :: w
    logical, intent(in) :: hover
    logical :: viewmode_activate_picking

    viewmode_activate_picking = .false.
    if (.not.hover) return

    if (w%viewmode < 0) then
       ! in window_forced mode, any click requires picking
       viewmode_activate_picking = any_mouse_clicked()
    elseif (w%viewmode == vm_navigate) then
       ! navigate -> when measuring, or on a double click (to clear the selection)
       viewmode_activate_picking = is_bind_event(BIND_NAV_MEASURE) .or.&
          is_bind_event(BIND_SELECT_MOLECULES_AND_DESELECT)
    elseif (w%viewmode == vm_select .or. w%viewmode == vm_movemol .or. w%viewmode == vm_moveatom) then
       ! select -> on any click, so the atom under the mouse is fresh
       ! move atoms -> on any click, to latch the grabbed atom/molecule
       viewmode_activate_picking = any_mouse_clicked()
    end if

  end function viewmode_activate_picking

  !> Process the mouse events in the view window, according to the
  !> different view modes. hover = whether the view is active and
  !> being hovered.
  module subroutine viewmode_process_events(w,hover)
    use interfaces_cimgui
    use scenes, only: scene
    use utils, only: translate, rotate, mult, invmult
    use tools_math, only: cross_cfloat, matinv_cfloat, axisangle2mat
    use keybindings, only: is_bind_event, is_bind_mousescroll, BIND_NAV_ROTATE,&
       BIND_NAV_ROTATE_PERP,&
       BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET, BIND_NAV_MEASURE,&
       BIND_CLOSE_FOCUSED_DIALOG, BIND_SELECT_MOLECULES_AND_DESELECT, BIND_SELECT_ATOMS,&
       BIND_SELECT_MOLECULES, BIND_MOVEMOL_TRANSLATE, BIND_MOVEMOL_ROTATE,&
       BIND_MOVEMOL_ROTATE_PERP, BIND_MOVEATOM_TRANSLATE
    use systems, only: nsys, sysc, sys, atlisttype_ncel_frac, lastchange_geometry
    use global, only: iunit_bohr
    use gui_main, only: io, ColorHighlightSelectScene
    class(window), intent(inout), target :: w
    logical, intent(in) :: hover

    type(ImVec2) :: texpos, mousepos, pmin, pmax
    real(c_float) :: ratio, pos3(3), vnew(3), vold(3), axis(3)
    real(c_float) :: mpos2(2), ang, xc(3), dist, comc(3)
    real*8 :: dxbohr(3)
    integer :: isys
    integer(c_int) :: col
    logical :: ok, dragged

    real(c_float), parameter :: mousesens_zoom0 = 0.15_c_float
    real(c_float), parameter :: mousesens_rot0 = 3._c_float
    real(c_float), parameter :: mousesens_vol0 = 0.05_c_float ! per-notch fractional cell-volume change
    real(c_float), parameter :: selrect_thr = 4._c_float ! click/drag threshold (pixels)

    ! first pass when opened, reset the state
    if (w%firstpass) then
       w%mposlast%x = 0._c_float
       w%mposlast%y = 0._c_float
       w%mpos0_r = 0._c_float
       w%mpos0_l = 0._c_float
       w%oldview = 0._c_float
       w%cpos0_l = 0._c_float
       w%mpos0_s = 0._c_float
       w%mpos0_m = 0._c_float
       w%ilock = ilock_no
       w%selrect_active = .false.
    end if

    ! window_forced view mode (pick atom): exits on a mouse click on the view
    ! (any button) or any key press anywhere; if an atom is under the mouse
    ! when clicked, save it as the pick result. Also exits if the window that
    ! commanded the mode is gone.
    if (w%viewmode == vm_pick_atom) then
       ! check the commanding window is still active
       ok = (w%vmdata%owner >= 1 .and. w%vmdata%owner <= nwin)
       if (ok) ok = win(w%vmdata%owner)%isinit .and. win(w%vmdata%owner)%isopen
       if (.not.ok) then
          w%vmdata%idx = 0
          w%viewmode = vm_navigate
          w%viewmode_transient = .false.
          return
       end if

       if (hover .and. any_mouse_clicked()) then
          ! pick the atom under the mouse, if any, and exit
          w%vmdata%idx = 0
          if (w%mousepos_idx(1) > 0) w%vmdata%idx = w%mousepos_idx
          w%viewmode = vm_navigate
          w%viewmode_transient = .false.
       elseif (.not.io%WantTextInput .and. any_key_pressed()) then
          ! cancelled; the result stays zero
          w%viewmode = vm_navigate
          w%viewmode_transient = .false.
          w%vmdata%idx = 0
       end if
       return
    end if

    ! only process if there is an associated system is viewed and scene is initialized
    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    if (.not.associated(w%sc)) return
    if (w%sc%isinit < 2) return

    ! release any interactive-dynamics grab whenever the view is not in the
    ! forced MD mode
    if (w%viewmode /= vm_mddrag) then
       if (w%ilock == ilock_mddrag .or. w%ilock == ilock_mdmovemol .or. w%ilock == ilock_mdrotmol) &
          w%ilock = ilock_no
       sysc(w%view_selected)%md%drag_iat = 0
       sysc(w%view_selected)%md%interacting = .false.
    end if

    ! process mode-specific events
    if (w%viewmode == vm_navigate .or. w%viewmode == vm_mddrag) then
       ! navigation and the forced interactive-dynamics mode share the camera
       ! controls below; vm_mddrag additionally grabs atoms/molecules first.
       isys = w%view_selected

       ! drop any rubber-band drag carried over from select mode (e.g. shift released mid-drag)
       w%selrect_active = .false.

       call igGetMousePos(mousepos)
       texpos = mousepos

       ! transform to the texture pos
       call w%mousepos_to_texpos(texpos)

       ! Interactive dynamics grabs (only in the forced MD
       ! mode). These run before the camera controls and take priority
       ! when the cursor is over an atom/molecule: left grabs an atom,
       ! right rigidly translates the molecule under the cursor,
       ! middle rigidly rotates it. When the cursor is on empty space
       ! the grab does not latch and the matching camera control below
       ! (rotate/translate/perp-rotate) runs instead.
       if (w%viewmode == vm_mddrag) then
          call md_atom_drag()
          call md_mol_move()
          call md_mol_rotate()
       end if

       ! Zoom. There are two behaviors: mouse scroll and hold key and
       ! translate mouse
       ratio = 0._c_float
       if (hover.and.(w%ilock == ilock_no .or. w%ilock == ilock_scroll).and. is_bind_event(BIND_NAV_ZOOM,.false.)) then
          if (is_bind_mousescroll(BIND_NAV_ZOOM)) then
             ! mouse scroll
             ratio = mousesens_zoom0 * io%MouseWheel
          else
             ! keys
             w%mpos0_s = mousepos%y
             w%ilock = ilock_scroll
          end if
       elseif (w%ilock == ilock_scroll) then
          if (is_bind_event(BIND_NAV_ZOOM,.true.)) then
             ! 10/a to make it adimensional
             ratio = mousesens_zoom0 * (w%mpos0_s-mousepos%y) * (10._c_float / w%FBOside)
             w%mpos0_s = mousepos%y
          else
             w%ilock = ilock_no
          end if
       end if
       if (ratio /= 0._c_float) then
          ratio = min(max(ratio,-0.99999_c_float),0.9999_c_float)
          call w%sc%cam_zoom(ratio)
          w%forcerender = .true.
       end if

       ! drag
       if (hover.and.is_bind_event(BIND_NAV_TRANSLATE,.false.).and.(w%ilock == ilock_no.or.w%ilock == ilock_right)) then
          ! drag on the scene-center depth plane; unprojecting at the far plane
          ! (z=1) is degenerate in perspective (point at infinity -> division by
          ! zero in unproject)
          vnew = w%sc%scenecenter
          call w%world_to_texpos(vnew)
          w%mpos0_r = (/texpos%x,texpos%y,vnew(3)/)

          ! save the current view matrix
          w%oldview = w%sc%view

          w%ilock = ilock_right
          w%mposlast = mousepos
       elseif (w%ilock == ilock_right) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (is_bind_event(BIND_NAV_TRANSLATE,.true.)) then
             if (mousepos%x /= w%mposlast%x .or. mousepos%y /= w%mposlast%y) then
                vnew = (/texpos%x,texpos%y,w%mpos0_r(3)/)
                call w%texpos_to_view(vnew)
                vold = w%mpos0_r
                call w%texpos_to_view(vold)

                xc = vold - vnew
                call invmult(xc,w%oldview)
                call w%sc%cam_move(xc)
                w%forcerender = .true.
             end if
          else
             w%ilock = ilock_no
          end if
       end if

       ! rotate
       if (hover .and. is_bind_event(BIND_NAV_ROTATE,.false.) .and. (w%ilock == ilock_no .or. w%ilock == ilock_left)) then
          w%mpos0_l = (/texpos%x, texpos%y, 0._c_float/)
          w%cpos0_l = w%mpos0_l
          call w%texpos_to_view(w%cpos0_l)
          w%ilock = ilock_left
       elseif (w%ilock == ilock_left) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (is_bind_event(BIND_NAV_ROTATE,.true.)) then
             if (texpos%x /= w%mpos0_l(1) .or. texpos%y /= w%mpos0_l(2)) then
                ! calculate the axis
                vnew = (/texpos%x,texpos%y,w%mpos0_l(3)/)
                call w%texpos_to_view(vnew)
                pos3 = (/0._c_float,0._c_float,1._c_float/)
                axis = cross_cfloat(pos3,vnew - w%cpos0_l)

                ! calculate the angle
                mpos2(1) = texpos%x - w%mpos0_l(1)
                mpos2(2) = texpos%y - w%mpos0_l(2)
                ang = 2._c_float * norm2(mpos2) * mousesens_rot0 / w%FBOside

                ! rotate the camera
                call w%sc%cam_rotate(axis,ang)

                ! save the new position and re-render
                w%mpos0_l = (/texpos%x, texpos%y, 0._c_float/)
                w%cpos0_l = w%mpos0_l
                call w%texpos_to_view(w%cpos0_l)
                w%forcerender = .true.
             end if
          else
             w%ilock = ilock_no
          end if
       end if

       ! rotate around perpendicular axis
       if (hover .and. is_bind_event(BIND_NAV_ROTATE_PERP,.false.) .and.&
          (w%ilock == ilock_no .or. w%ilock == ilock_middle)) then
          w%mpos0_m = w%sc%scenecenter
          call w%world_to_texpos(w%mpos0_m)
          vnew = (/texpos%x,texpos%y,0._c_float/)
          vnew = vnew - w%mpos0_m
          dist = norm2(vnew)
          if (dist > 0._c_float) then
             w%cpos0_m = vnew / dist
          else
             w%cpos0_m = (/0._c_float, -1._c_float, 0._c_float/)
          end if
          w%ilock = ilock_middle
       elseif (w%ilock == ilock_middle) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (is_bind_event(BIND_NAV_ROTATE_PERP,.true.)) then
             vnew = (/texpos%x,texpos%y,0._c_float/)
             vnew = vnew - w%mpos0_m
             dist = norm2(vnew)
             if (dist > 0._c_float) then
                vnew = vnew / dist
                xc = cross_cfloat(w%cpos0_m,vnew)
                ang = atan2(xc(3),dot_product(w%cpos0_m,vnew))
                axis = (/0._c_float,0._c_float,1._c_float/)
                call w%sc%cam_rotate(axis,ang)
                w%forcerender = .true.
                w%cpos0_m = vnew
             end if
          else
             w%ilock = ilock_no
          end if
       end if

       ! reset the view
       if (hover .and. is_bind_event(BIND_NAV_RESET,.false.)) then
          call w%sc%reset()
          w%forcerender = .true.
       end if

       ! atom selection
       if (hover .and. is_bind_event(BIND_NAV_MEASURE)) then
          call w%sc%select_atom(w%mousepos_idx)
          w%forcerender = .true.
       end if
       if (hover .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) then
          call w%sc%select_atom((/0,0,0,0,0/))
          w%forcerender = .true.
       end if

       ! double click on empty space clears the selection
       if (hover .and. igIsMouseDoubleClicked(ImGuiMouseButton_Left) .and. w%mousepos_idx(1) == 0) then
          call sysc(w%view_selected)%highlight_clear(.false.)
          w%forcerender = .true.
       end if
    elseif (w%viewmode == vm_select) then
       ! select mode
       isys = w%view_selected
       call igGetMousePos(mousepos)

       ! select the whole fragment under the mouse with either a right click or
       ! a left double click; a single left click starts a potential rubber band
       if (hover) then
          if (is_bind_event(BIND_SELECT_MOLECULES_AND_DESELECT)) then
             if (w%mousepos_idx(1) > 0) then
                ! the first click of the pair already toggled this single atom on
                ! its release; undo it before toggling the fragment
                call toggle_cellatoms((/w%mousepos_idx(1)/))
                call toggle_fragment(w%mousepos_idx(1))
             else
                ! double click on empty space clears the selection
                call sysc(isys)%highlight_clear(.false.)
             end if
             w%selrect_active = .false.
             w%forcerender = .true.
          elseif (is_bind_event(BIND_SELECT_MOLECULES) .and. w%mousepos_idx(1) > 0) then
             call toggle_fragment(w%mousepos_idx(1))
             w%selrect_active = .false.
             w%forcerender = .true.
          elseif (is_bind_event(BIND_SELECT_ATOMS)) then
             ! left click: start a potential rubber-band drag
             w%selrect_active = .true.
             w%selrect_p0 = mousepos
          end if
       end if

       ! left drag/release: rubber-band rectangle or single-atom toggle
       if (w%selrect_active) then
          dragged = abs(mousepos%x-w%selrect_p0%x) > selrect_thr .or.&
             abs(mousepos%y-w%selrect_p0%y) > selrect_thr
          if (is_bind_event(BIND_SELECT_ATOMS,.true.)) then
             ! draw the rubber band once dragged beyond the click/drag threshold
             if (dragged) then
                pmin%x = min(w%selrect_p0%x,mousepos%x)
                pmin%y = min(w%selrect_p0%y,mousepos%y)
                pmax%x = max(w%selrect_p0%x,mousepos%x)
                pmax%y = max(w%selrect_p0%y,mousepos%y)
                ! translucent fill plus a slightly darker border
                col = igGetColorU32_Vec4(ImVec4(ColorHighlightSelectScene(1),&
                   ColorHighlightSelectScene(2),ColorHighlightSelectScene(3),0.3_c_float))
                call ImDrawList_AddRectFilled(igGetWindowDrawList(),pmin,pmax,col,0._c_float,0_c_int)
                col = igGetColorU32_Vec4(ImVec4(0.6_c_float*ColorHighlightSelectScene(1),&
                   0.6_c_float*ColorHighlightSelectScene(2),0.6_c_float*ColorHighlightSelectScene(3),0.8_c_float))
                call ImDrawList_AddRect(igGetWindowDrawList(),pmin,pmax,col,0._c_float,0_c_int,1.5_c_float)
             end if
          else
             ! mouse released: a drag selects a rectangle, a static click toggles one atom
             if (dragged) then
                call select_in_rect(w%selrect_p0,mousepos)
             elseif (w%mousepos_idx(1) > 0) then
                call toggle_cellatoms((/w%mousepos_idx(1)/))
             end if
             w%selrect_active = .false.
             w%forcerender = .true.
          end if
       end if
    elseif (w%viewmode == vm_movemol) then
       ! move-molecules mode: rigidly drag the molecule/atom under the cursor,
       ! preserving the bonding
       isys = w%view_selected
       call igGetMousePos(mousepos)
       texpos = mousepos
       call w%mousepos_to_texpos(texpos)

       ! scroll: resize the cell (crystal) or zoom the camera (molecule)
       call moveobj_scroll()

       ! translate (right mouse): whole molecule if the fragment is discrete,
       ! otherwise just the single atom; the grabbed atom stays under the cursor
       if (hover.and.is_bind_event(BIND_MOVEMOL_TRANSLATE,.false.).and.&
          (w%ilock == ilock_no.or.w%ilock == ilock_right)) then
          call moveobj_latch()
          if (w%moveobj_icel > 0) then
             ! drag on the scene-center depth plane (see note above): the far
             ! plane (z=1) is at infinity in perspective and crashes unproject
             vnew = w%sc%scenecenter
             call w%world_to_texpos(vnew)
             w%mpos0_r = (/texpos%x,texpos%y,vnew(3)/)
             w%ilock = ilock_right
             w%mposlast = mousepos
          end if
       elseif (w%ilock == ilock_right) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (w%moveobj_icel > 0 .and. is_bind_event(BIND_MOVEMOL_TRANSLATE,.true.)) then
             if (mousepos%x /= w%mposlast%x .or. mousepos%y /= w%mposlast%y) then
                ! world (bohr) displacement matching the cursor motion
                vnew = (/texpos%x,texpos%y,w%mpos0_r(3)/)
                call w%texpos_to_view(vnew)
                vold = w%mpos0_r
                call w%texpos_to_view(vold)
                xc = vnew - vold
                call invmult(xc,w%sc%view,notrans=.true.)  ! eye -> tworld
                call invmult(xc,w%sc%world,notrans=.true.) ! tworld -> world (bohr)
                dxbohr = real(xc,8)
                if (w%moveobj_isdiscrete) then
                   call sys(isys)%c%move_molecule(w%moveobj_imol,dxbohr,iunit_bohr,&
                      .true.,copybonding=.true.)
                else
                   call sys(isys)%c%move_atom(w%moveobj_icel,dxbohr,iunit_bohr,&
                      .false.,.true.,copybonding=.true.)
                end if
                sysc(isys)%sc%nextbuildlists_fixcam = .true.
                call sysc(isys)%post_event(lastchange_geometry)
                w%forcerender = .true.
                w%mpos0_r = (/texpos%x,texpos%y,w%mpos0_r(3)/)
                w%mposlast = mousepos
             end if
          else
             w%ilock = ilock_no
          end if
       end if

       ! rotate the molecule about its COM (left mouse), discrete fragments only
       if (hover.and.is_bind_event(BIND_MOVEMOL_ROTATE,.false.).and.&
          (w%ilock == ilock_no.or.w%ilock == ilock_left)) then
          call moveobj_latch()
          if (w%moveobj_icel > 0 .and. w%moveobj_isdiscrete) then
             w%mpos0_l = (/texpos%x,texpos%y,0._c_float/)
             w%cpos0_l = w%mpos0_l
             call w%texpos_to_view(w%cpos0_l)
             w%ilock = ilock_left
          end if
       elseif (w%ilock == ilock_left) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (w%moveobj_icel > 0 .and. w%moveobj_isdiscrete .and.&
             is_bind_event(BIND_MOVEMOL_ROTATE,.true.)) then
             if (texpos%x /= w%mpos0_l(1) .or. texpos%y /= w%mpos0_l(2)) then
                ! arcball axis (eye) + angle, same math as scene rotation
                vnew = (/texpos%x,texpos%y,w%mpos0_l(3)/)
                call w%texpos_to_view(vnew)
                pos3 = (/0._c_float,0._c_float,1._c_float/)
                axis = cross_cfloat(pos3,vnew - w%cpos0_l)
                mpos2(1) = texpos%x - w%mpos0_l(1)
                mpos2(2) = texpos%y - w%mpos0_l(2)
                ang = 2._c_float * norm2(mpos2) * mousesens_rot0 / w%FBOside
                call movemol_rotate_molecule(axis,ang)
                w%forcerender = .true.
                w%mpos0_l = (/texpos%x,texpos%y,0._c_float/)
                w%cpos0_l = w%mpos0_l
                call w%texpos_to_view(w%cpos0_l)
             end if
          else
             w%ilock = ilock_no
          end if
       end if

       ! rotate the molecule about the screen-perpendicular axis (middle mouse),
       ! discrete fragments only
       if (hover.and.is_bind_event(BIND_MOVEMOL_ROTATE_PERP,.false.).and.&
          (w%ilock == ilock_no.or.w%ilock == ilock_middle)) then
          call moveobj_latch()
          if (w%moveobj_icel > 0 .and. w%moveobj_isdiscrete) then
             ! project the molecule center of mass to screen
             comc = sys(isys)%c%mol(w%moveobj_imol)%cmass()
             call mult(w%mpos0_m,w%sc%world,comc)
             call w%world_to_texpos(w%mpos0_m)
             vnew = (/texpos%x,texpos%y,0._c_float/)
             vnew = vnew - w%mpos0_m
             dist = norm2(vnew)
             if (dist > 0._c_float) then
                w%cpos0_m = vnew / dist
             else
                w%cpos0_m = (/0._c_float,-1._c_float,0._c_float/)
             end if
             w%ilock = ilock_middle
          end if
       elseif (w%ilock == ilock_middle) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (w%moveobj_icel > 0 .and. w%moveobj_isdiscrete .and.&
             is_bind_event(BIND_MOVEMOL_ROTATE_PERP,.true.)) then
             vnew = (/texpos%x,texpos%y,0._c_float/)
             vnew = vnew - w%mpos0_m
             dist = norm2(vnew)
             if (dist > 0._c_float) then
                vnew = vnew / dist
                xc = cross_cfloat(w%cpos0_m,vnew)
                ang = atan2(xc(3),dot_product(w%cpos0_m,vnew))
                axis = (/0._c_float,0._c_float,1._c_float/)
                call movemol_rotate_molecule(axis,ang)
                w%forcerender = .true.
                w%cpos0_m = vnew
             end if
          else
             w%ilock = ilock_no
          end if
       end if
    elseif (w%viewmode == vm_moveatom) then
       ! move-atoms mode: left- or right-drag translates the single atom under
       ! the cursor (never the whole molecule); scroll resizes the cell (crystals)
       isys = w%view_selected
       call igGetMousePos(mousepos)
       texpos = mousepos
       call w%mousepos_to_texpos(texpos)

       ! scroll: resize the cell (crystal) or zoom the camera (molecule)
       call moveobj_scroll()

       ! translate the single atom under the cursor (left mouse)
       call moveatom_translate(BIND_MOVEATOM_TRANSLATE,ilock_left,w%mpos0_l)
    end if

    ! if this is a transient view mode, reset to default (navigation)
    if (w%viewmode_transient) then
       w%viewmode = vm_navigate
       w%viewmode_transient = .false.
    end if

  contains
    ! whether any (non-modifier) key was pressed this frame
    function any_key_pressed()
      use keybindings, only: is_mod_key
      logical :: any_key_pressed
      integer :: k

      any_key_pressed = .false.
      do k = ImGuiKey_NamedKey_BEGIN, ImGuiKey_NamedKey_END-1
         if (is_mod_key(k)) cycle
         if (igIsKeyPressed(k,.false._c_bool)) then
            any_key_pressed = .true.
            return
         end if
      end do
    end function any_key_pressed

    ! mouse scroll while in a move mode: change the cell volume isotropically
    ! for a crystal, otherwise fall back to the camera zoom
    subroutine moveobj_scroll()
      real(c_float) :: rr

      if (.not.hover) return
      if (.not.(w%ilock == ilock_no .or. w%ilock == ilock_scroll)) return
      if (io%MouseWheel == 0._c_float) return
      if (sys(isys)%c%ismolecule) then
         rr = min(max(mousesens_zoom0*io%MouseWheel,-0.99999_c_float),0.9999_c_float)
         call w%sc%cam_zoom(rr)
         w%forcerender = .true.
      else
         ! rescale the cell volume but keep the current bonds (copybonding)
         call sys(isys)%c%move_cell(0,1d0 + real(mousesens_vol0*io%MouseWheel,8),iunit_bohr,&
            .false.,.true.,copybonding=.true.)
         sysc(isys)%sc%nextbuildlists_fixcam = .true.
         call sysc(isys)%post_event(lastchange_geometry)
         w%forcerender = .true.
      end if
    end subroutine moveobj_scroll

    ! translate the single latched atom by dragging with the given mouse bind;
    ! ilockval and anchor are the per-button lock state and drag anchor. Always
    ! moves just the one atom (used by the move-atoms mode)
    subroutine moveatom_translate(bindid,ilockval,anchor)
      integer, intent(in) :: bindid, ilockval
      real(c_float), intent(inout) :: anchor(3)

      if (hover.and.is_bind_event(bindid,.false.).and.&
         (w%ilock == ilock_no.or.w%ilock == ilockval)) then
         call moveobj_latch()
         if (w%moveobj_icel > 0) then
            ! drag on the scene-center depth plane (the far plane is at infinity
            ! in perspective and crashes unproject)
            vnew = w%sc%scenecenter
            call w%world_to_texpos(vnew)
            anchor = (/texpos%x,texpos%y,vnew(3)/)
            w%ilock = ilockval
            w%mposlast = mousepos
         end if
      elseif (w%ilock == ilockval) then
         call igSetMouseCursor(ImGuiMouseCursor_Hand)
         if (w%moveobj_icel > 0 .and. is_bind_event(bindid,.true.)) then
            if (mousepos%x /= w%mposlast%x .or. mousepos%y /= w%mposlast%y) then
               ! world (bohr) displacement matching the cursor motion
               vnew = (/texpos%x,texpos%y,anchor(3)/)
               call w%texpos_to_view(vnew)
               vold = anchor
               call w%texpos_to_view(vold)
               xc = vnew - vold
               call invmult(xc,w%sc%view,notrans=.true.)  ! eye -> tworld
               call invmult(xc,w%sc%world,notrans=.true.) ! tworld -> world (bohr)
               dxbohr = real(xc,8)
               call sys(isys)%c%move_atom(w%moveobj_icel,dxbohr,iunit_bohr,&
                  .false.,.true.,copybonding=.true.)
               sysc(isys)%sc%nextbuildlists_fixcam = .true.
               call sysc(isys)%post_event(lastchange_geometry)
               w%forcerender = .true.
               anchor = (/texpos%x,texpos%y,anchor(3)/)
               w%mposlast = mousepos
            end if
         else
            w%ilock = ilock_no
         end if
      end if
    end subroutine moveatom_translate

    ! Interactive dynamics (vm_mddrag): grab the atom under the cursor and
    ! drag it to follow the mouse (left button), by setting the MD drag
    ! target that the integrator clamps each step. Latches on the grabbed atom's
    ! depth plane. Does nothing (leaving the camera rotation to run) if the
    ! cursor is not over an atom.
    subroutine md_atom_drag()
      if (hover.and.is_bind_event(BIND_NAV_ROTATE,.false.).and.w%ilock == ilock_no.and.&
         w%mousepos_idx(1) > 0 .and. sysc(isys)%md%ready) then
         sysc(isys)%md%drag_iat = w%mousepos_idx(1)
         sysc(isys)%md%drag_target = sys(isys)%c%atcel(w%mousepos_idx(1))%r
         vnew = real(sys(isys)%c%atcel(w%mousepos_idx(1))%r,c_float)
         call w%world_to_texpos(vnew)
         w%mpos0_r = (/texpos%x,texpos%y,vnew(3)/)
         w%ilock = ilock_mddrag
         w%mposlast = mousepos
      elseif (w%ilock == ilock_mddrag) then
         call igSetMouseCursor(ImGuiMouseCursor_Hand)
         if (is_bind_event(BIND_NAV_ROTATE,.true.)) then
            if (mousepos%x /= w%mposlast%x .or. mousepos%y /= w%mposlast%y) then
               call drag_delta_world(dxbohr)
               sysc(isys)%md%drag_target = sysc(isys)%md%drag_target + dxbohr
               w%mpos0_r = (/texpos%x,texpos%y,w%mpos0_r(3)/)
               w%mposlast = mousepos
               w%forcerender = .true.
            end if
         else
            w%ilock = ilock_no
            sysc(isys)%md%drag_iat = 0
         end if
      end if
    end subroutine md_atom_drag

    ! interactive dynamics (vm_mddrag): rigidly translate the molecule under the
    ! cursor (right button) so it follows the mouse, by shifting its atoms in the
    ! MD position buffer. Does nothing (leaving the camera pan to run) if the
    ! cursor is not over an atom.
    subroutine md_mol_move()
      if (hover.and.is_bind_event(BIND_NAV_TRANSLATE,.false.).and.w%ilock == ilock_no.and.&
         w%mousepos_idx(1) > 0 .and. sysc(isys)%md%ready) then
         call moveobj_latch()
         if (w%moveobj_icel > 0) then
            ! drag on the scene-center depth plane (far plane -> unproject crash)
            vnew = w%sc%scenecenter
            call w%world_to_texpos(vnew)
            w%mpos0_r = (/texpos%x,texpos%y,vnew(3)/)
            w%ilock = ilock_mdmovemol
            w%mposlast = mousepos
            sysc(isys)%md%interacting = .true.
         end if
      elseif (w%ilock == ilock_mdmovemol) then
         call igSetMouseCursor(ImGuiMouseCursor_Hand)
         if (w%moveobj_icel > 0 .and. is_bind_event(BIND_NAV_TRANSLATE,.true.)) then
            if (mousepos%x /= w%mposlast%x .or. mousepos%y /= w%mposlast%y) then
               call drag_delta_world(dxbohr)
               call md_move_fragment(dxbohr)
               w%mpos0_r = (/texpos%x,texpos%y,w%mpos0_r(3)/)
               w%mposlast = mousepos
               w%forcerender = .true.
            end if
         else
            w%ilock = ilock_no
            sysc(isys)%md%interacting = .false.
         end if
      end if
    end subroutine md_mol_move

    ! interactive dynamics (vm_mddrag): rigidly rotate the molecule under the
    ! cursor about its center of mass (middle button), arcball-style, by
    ! rotating its atoms in the MD position buffer. Does nothing (leaving the
    ! camera perpendicular rotation to run) if the cursor is not over an atom.
    subroutine md_mol_rotate()
      real(c_float) :: axisw(3), lax
      real*8 :: rmat(3,3)

      if (hover.and.is_bind_event(BIND_NAV_ROTATE_PERP,.false.).and.w%ilock == ilock_no.and.&
         w%mousepos_idx(1) > 0 .and. sysc(isys)%md%ready) then
         call moveobj_latch()
         if (w%moveobj_icel > 0) then
            w%mpos0_m = (/texpos%x,texpos%y,0._c_float/)
            w%cpos0_m = w%mpos0_m
            call w%texpos_to_view(w%cpos0_m)
            w%ilock = ilock_mdrotmol
            sysc(isys)%md%interacting = .true.
         end if
      elseif (w%ilock == ilock_mdrotmol) then
         call igSetMouseCursor(ImGuiMouseCursor_Hand)
         if (w%moveobj_icel > 0 .and. is_bind_event(BIND_NAV_ROTATE_PERP,.true.)) then
            if (texpos%x /= w%mpos0_m(1) .or. texpos%y /= w%mpos0_m(2)) then
               ! arcball axis (eye) + angle, same math as the camera rotation
               vnew = (/texpos%x,texpos%y,w%mpos0_m(3)/)
               call w%texpos_to_view(vnew)
               pos3 = (/0._c_float,0._c_float,1._c_float/)
               axis = cross_cfloat(pos3,vnew - w%cpos0_m)
               mpos2(1) = texpos%x - w%mpos0_m(1)
               mpos2(2) = texpos%y - w%mpos0_m(2)
               ang = 2._c_float * norm2(mpos2) * mousesens_rot0 / w%FBOside
               ! axis from eye to world (bohr), then rotate the fragment
               lax = norm2(axis)
               if (lax > 1e-10_c_float) then
                  axisw = axis / lax
                  call invmult(axisw,w%sc%world,notrans=.true.)
                  if (norm2(axisw) > 1e-10_c_float) then
                     rmat = axisangle2mat(real(axisw,8),real(ang,8))
                     call md_rotate_fragment(rmat)
                     w%forcerender = .true.
                  end if
               end if
               w%mpos0_m = (/texpos%x,texpos%y,0._c_float/)
               w%cpos0_m = w%mpos0_m
               call w%texpos_to_view(w%cpos0_m)
            end if
         else
            w%ilock = ilock_no
            sysc(isys)%md%interacting = .false.
         end if
      end if
    end subroutine md_mol_rotate

    ! world displacement matching the cursor motion from the drag anchor
    ! w%mpos0_r on its depth plane (used by the MD atom/molecule translation drags)
    subroutine drag_delta_world(dbohr)
      real*8, intent(out) :: dbohr(3)
      real(c_float) :: va(3), vb(3), dc(3)

      va = (/texpos%x,texpos%y,w%mpos0_r(3)/)
      call w%texpos_to_view(va)
      vb = w%mpos0_r
      call w%texpos_to_view(vb)
      dc = va - vb
      call invmult(dc,w%sc%view,notrans=.true.)  ! eye -> tworld
      call invmult(dc,w%sc%world,notrans=.true.) ! tworld -> world (bohr)
      dbohr = real(dc,8)
    end subroutine drag_delta_world

    ! shift the latched fragment's atoms in the MD position buffer by dx (bohr):
    ! the whole molecule if it is discrete, otherwise just the single latched
    ! atom (matching the move-molecules mode); keep c synced
    subroutine md_move_fragment(dx)
      real*8, intent(in) :: dx(3)
      integer :: imol, k, iat

      imol = w%moveobj_imol
      if (w%moveobj_isdiscrete .and. imol >= 1 .and. imol <= sys(isys)%c%nmol) then
         do k = 1, sys(isys)%c%mol(imol)%nat
            iat = sys(isys)%c%mol(imol)%at(k)%cidx
            sysc(isys)%md%r(:,iat) = sysc(isys)%md%r(:,iat) + dx
         end do
      elseif (w%moveobj_icel > 0) then
         sysc(isys)%md%r(:,w%moveobj_icel) = sysc(isys)%md%r(:,w%moveobj_icel) + dx
      end if
      call sys(isys)%c%update_positions(sysc(isys)%md%r)
    end subroutine md_move_fragment

    ! rotate the latched molecule's atoms in the MD position buffer about their
    ! mass-weighted center by rmat (world/bohr); discrete molecules only (as in
    ! the move-molecules mode); keep c synced
    subroutine md_rotate_fragment(rmat)
      real*8, intent(in) :: rmat(3,3)
      integer :: imol, k, iat
      real*8 :: com(3), mtot

      imol = w%moveobj_imol
      if (.not.w%moveobj_isdiscrete .or. imol < 1 .or. imol > sys(isys)%c%nmol) return
      com = 0d0
      mtot = 0d0
      do k = 1, sys(isys)%c%mol(imol)%nat
         iat = sys(isys)%c%mol(imol)%at(k)%cidx
         com = com + sysc(isys)%md%mass(iat) * sysc(isys)%md%r(:,iat)
         mtot = mtot + sysc(isys)%md%mass(iat)
      end do
      if (mtot <= 0d0) return
      com = com / mtot
      do k = 1, sys(isys)%c%mol(imol)%nat
         iat = sys(isys)%c%mol(imol)%at(k)%cidx
         sysc(isys)%md%r(:,iat) = com + matmul(rmat, sysc(isys)%md%r(:,iat) - com)
      end do
      call sys(isys)%c%update_positions(sysc(isys)%md%r)
    end subroutine md_rotate_fragment

    ! latch the cell atom and its molecular fragment under the cursor at the
    ! start of a move drag
    subroutine moveobj_latch()
      integer :: jmol

      w%moveobj_icel = 0
      w%moveobj_imol = 0
      w%moveobj_isdiscrete = .false.
      if (w%mousepos_idx(1) <= 0) return
      w%moveobj_icel = w%mousepos_idx(1)
      if (allocated(sys(isys)%c%idatcelmol)) then
         jmol = sys(isys)%c%idatcelmol(1,w%moveobj_icel)
         w%moveobj_imol = jmol
         if (jmol >= 1 .and. jmol <= sys(isys)%c%nmol) &
            w%moveobj_isdiscrete = sys(isys)%c%mol(jmol)%discrete
      end if
    end subroutine moveobj_latch

    ! rotate the latched molecule rigidly about its center of mass,
    ! preserving bonding. axis0 is in eye/view coordinates (as
    ! produced by the navigation arcball math); ang0 is the rotation
    ! angle (radians).
    subroutine movemol_rotate_molecule(axis0,ang0)
      use tools_math, only: euler2mat, axisangle2mat
      real(c_float), intent(in) :: axis0(3)
      real(c_float), intent(in) :: ang0

      real(c_float) :: axisw(3), lax
      real*8 :: rinc(3,3)
      integer :: imol

      imol = w%moveobj_imol
      if (imol < 1) return

      ! axis from eye/view to world (bohr), matching scene_cam_rotate
      lax = norm2(axis0)
      if (lax <= 1e-10_c_float) return
      axisw = axis0 / lax
      call invmult(axisw,w%sc%world,notrans=.true.)
      if (norm2(axisw) <= 1e-10_c_float) return

      ! incremental rotation (Rodrigues) about the world-space axis, composed
      ! with the molecule's current standard-frame orientation
      rinc = axisangle2mat(real(axisw,8),real(ang0,8))
      if (.not.sys(isys)%c%mol(imol)%axes_computed) call sys(isys)%c%mol(imol)%compute_std()
      rinc = matmul(rinc,euler2mat(sys(isys)%c%mol(imol)%euler_std))
      call sys(isys)%c%rotate_molecule(imol,rmat=rinc,copybonding=.true.)
      sysc(isys)%sc%nextbuildlists_fixcam = .true.
      call sysc(isys)%post_event(lastchange_geometry)

    end subroutine movemol_rotate_molecule

    ! whether cell atom icel is currently in the persistent selection
    function cellatom_selected(icel)
      logical :: cellatom_selected
      integer, intent(in) :: icel

      cellatom_selected = .false.
      if (allocated(sysc(isys)%highlight_rgba)) then
         if (icel >= 1 .and. icel <= size(sysc(isys)%highlight_rgba,2)) &
            cellatom_selected = all(sysc(isys)%highlight_rgba(:,icel) >= 0._c_float)
      end if
    end function cellatom_selected

    ! toggle a list of cell atoms: if all are already selected, clear them;
    ! otherwise select all of them
    subroutine toggle_cellatoms(idlist)
      integer, intent(in) :: idlist(:)

      integer :: k
      logical :: allsel

      if (size(idlist) == 0) return
      allsel = .true.
      do k = 1, size(idlist)
         if (.not.cellatom_selected(idlist(k))) then
            allsel = .false.
            exit
         end if
      end do
      if (allsel) then
         call sysc(isys)%highlight_clear(.false.,idlist,atlisttype_ncel_frac)
      else
         call sysc(isys)%highlight_atoms(.false.,idlist,atlisttype_ncel_frac,&
            spread(ColorHighlightSelectScene,2,size(idlist)))
      end if
    end subroutine toggle_cellatoms

    ! toggle the whole fragment (molecule) that cell atom icel belongs to
    subroutine toggle_fragment(icel)
      integer, intent(in) :: icel

      integer :: jmol, k
      integer, allocatable :: flist(:)

      jmol = sys(isys)%c%idatcelmol(1,icel)
      allocate(flist(sys(isys)%c%mol(jmol)%nat))
      do k = 1, sys(isys)%c%mol(jmol)%nat
         flist(k) = sys(isys)%c%mol(jmol)%at(k)%cidx
      end do
      call toggle_cellatoms(flist)
    end subroutine toggle_fragment

    ! Add to the selection every atom whose projected center falls inside the
    ! screen-space rectangle delimited by p0 and p1.
    subroutine select_in_rect(p0,p1)
      type(ImVec2), intent(in) :: p0, p1

      type(ImVec2) :: t0, t1
      integer :: isph, jcel, n
      real(c_float) :: xmin, xmax, ymin, ymax, xp(3), xt(3)
      integer, allocatable :: lst(:)
      logical, allocatable :: seen(:)

      if (.not.associated(w%sc)) return
      if (w%sc%isinit < 2) return
      if (sys(isys)%c%ncel <= 0) return
      if (w%sc%obj%nsph <= 0) return

      ! rectangle corners: mouse position -> texture position
      t0 = p0
      t1 = p1
      call w%mousepos_to_texpos(t0)
      call w%mousepos_to_texpos(t1)
      xmin = min(t0%x,t1%x)
      xmax = max(t0%x,t1%x)
      ymin = min(t0%y,t1%y)
      ymax = max(t0%y,t1%y)

      ! project every drawn sphere center and collect the unique cell atoms whose
      ! centers fall inside the rectangle
      allocate(lst(w%sc%obj%nsph),seen(sys(isys)%c%ncel))
      seen = .false.
      n = 0
      do isph = 1, w%sc%obj%nsph
         jcel = w%sc%obj%sph(isph)%idx(1)
         if (jcel < 1 .or. jcel > sys(isys)%c%ncel) cycle
         if (seen(jcel)) cycle

         ! world coordinates -> texture position (x,y) + depth (z)
         xt = w%sc%obj%sph(isph)%x
         call mult(xp,w%sc%world,xt)
         call w%world_to_texpos(xp)
         if (xp(3) < 0._c_float .or. xp(3) > 1._c_float) cycle ! outside the view frustum
         if (xp(1) < xmin .or. xp(1) > xmax) cycle
         if (xp(2) < ymin .or. xp(2) > ymax) cycle

         n = n + 1
         lst(n) = jcel
         seen(jcel) = .true.
      end do

      if (n > 0) &
         call sysc(isys)%highlight_atoms(.false.,lst(1:n),atlisttype_ncel_frac,&
            spread(ColorHighlightSelectScene,2,n))
    end subroutine select_in_rect
  end subroutine viewmode_process_events

  !> Whether any mouse button was clicked this frame
  function any_mouse_clicked()
    logical :: any_mouse_clicked

    any_mouse_clicked = igIsMouseClicked(ImGuiMouseButton_Left,.false._c_bool) .or.&
       igIsMouseClicked(ImGuiMouseButton_Right,.false._c_bool) .or.&
       igIsMouseClicked(ImGuiMouseButton_Middle,.false._c_bool)

  end function any_mouse_clicked

  module subroutine draw_selection_tooltip(w,idx)
    use interfaces_cimgui
    use utils, only: iw_text
    use systems, only: sys
    use gui_main, only: fontsize, ColorMeasureSelect
    use tools_io, only: string
    use param, only: bohrtoa, pi
    class(window), intent(inout), target :: w
    integer(c_int), intent(in) :: idx(5)

    integer :: nmsel
    integer :: msel(5,4)
    integer :: idx1(4), idx2(4), idx3(4), idx4(4)
    real*8 :: x0(3), d, d1, d2, ang

    if (.not.associated(w%sc)) return

    ! check if the tooltip is needed
    nmsel = w%sc%nmsel
    if (nmsel == 0) return
    msel = w%sc%msel
    if (nmsel == 1 .and. (idx(1) == 0 .or. idx(5) == msel(5,1))) return

    ! start tooltip and header
    call igBeginTooltip()
    call igPushTextWrapPos(60._c_float * fontsize%x)
    call iw_text("Distance (d), angle (α), dihedral (φ)")

    ! distance 1-2
    idx1 = msel(1:4,1)
    if (nmsel == 1) then
       idx2 = idx(1:4)
       if (idx(1) == 0) goto 999
    else
       idx2 = msel(1:4,2)
    end if
    x0 = sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4)
    x0 = x0 - (sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4))
    x0 = sys(w%view_selected)%c%x2c(x0)
    d = norm2(x0)*bohrtoa
    if (abs(d) > 1d-14) then
       call iw_text("d(")
       call iw_text("1",rgb=ColorMeasureSelect(1:3,1),sameline_nospace=.true.)
       if (nmsel > 1) then
          call iw_text("2",rgb=ColorMeasureSelect(1:3,2),sameline_nospace=.true.)
       else
          call iw_text("*",sameline_nospace=.true.)
       end if
       call iw_text(")=" // string(d,'f',decimal=4) // " Å",sameline_nospace=.true.)
    end if

    ! distance and angle with atom 3
    if (nmsel > 1) then
       ! distance 2-3
       idx1 = msel(1:4,2)
       if (nmsel == 2) then
          idx2 = idx(1:4)
          if (idx(1) == 0) goto 999
          if (idx(5) == msel(5,2) .or. idx(5) == msel(5,1)) goto 999
       else
          idx2 = msel(1:4,3)
       end if
       x0 = sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4)
       x0 = x0 - (sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4))
       x0 = sys(w%view_selected)%c%x2c(x0)
       d = norm2(x0)*bohrtoa
       if (d > 1d-14) then
          call iw_text("d(")
          call iw_text("2",rgb=ColorMeasureSelect(1:3,2),sameline_nospace=.true.)
          if (nmsel > 2) then
             call iw_text("3",rgb=ColorMeasureSelect(1:3,3),sameline_nospace=.true.)
          else
             call iw_text("*",sameline_nospace=.true.)
          end if
          call iw_text(")=" // string(d,'f',decimal=4) // " Å",sameline_nospace=.true.)
       end if

       ! angle 1-2-3
       idx3 = msel(1:4,1)
       d1 = sys(w%view_selected)%c%distance(sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4),&
          sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       d2 = sys(w%view_selected)%c%distance(sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4),&
          sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       if (d1 > 1d-14 .and. d2 > 1d-14) then
          ang = sys(w%view_selected)%c%angle(&
             sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4),&
             sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4),&
             sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4)) * 180d0 / pi
          call iw_text(", α(",sameline_nospace=.true.)
          call iw_text("1",rgb=ColorMeasureSelect(1:3,1),sameline_nospace=.true.)
          call iw_text("2",rgb=ColorMeasureSelect(1:3,2),sameline_nospace=.true.)
          if (nmsel > 2) then
             call iw_text("3",rgb=ColorMeasureSelect(1:3,3),sameline_nospace=.true.)
          else
             call iw_text("*",sameline_nospace=.true.)
          end if
          call iw_text(")=" // string(ang,'f',decimal=2) // "°",sameline_nospace=.true.)
       end if
    end if

    ! distance, angle, dihedral
    if (nmsel > 2) then
       ! distance 3-4
       idx1 = msel(1:4,3)
       if (nmsel == 3) then
          idx2 = idx(1:4)
          if (idx(1) == 0) goto 999
          if (idx(5) == msel(5,3) .or. idx(5) == msel(5,2) .or. idx(5) == msel(5,1)) goto 999
       else
          idx2 = msel(1:4,4)
       end if
       x0 = sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4)
       x0 = x0 - (sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4))
       x0 = sys(w%view_selected)%c%x2c(x0)
       d = norm2(x0)*bohrtoa
       if (d > 1d-14) then
          call iw_text("d(")
          call iw_text("3",rgb=ColorMeasureSelect(1:3,3),sameline_nospace=.true.)
          if (nmsel > 3) then
             call iw_text("4",rgb=ColorMeasureSelect(1:3,4),sameline_nospace=.true.)
          else
             call iw_text("*",sameline_nospace=.true.)
          end if
          call iw_text(")=" // string(d,'f',decimal=4) // " Å",sameline_nospace=.true.)
       end if

       ! angle 2-3-4
       idx3 = msel(1:4,2)
       d1 = sys(w%view_selected)%c%distance(sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4),&
          sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       d2 = sys(w%view_selected)%c%distance(sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4),&
          sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       if (d1 > 1d-14 .and. d2 > 1d-14) then
          ang = sys(w%view_selected)%c%angle(&
             sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4),&
             sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4),&
             sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4)) * 180d0 / pi
          call iw_text(", α(",sameline_nospace=.true.)
          call iw_text("2",rgb=ColorMeasureSelect(1:3,2),sameline_nospace=.true.)
          call iw_text("3",rgb=ColorMeasureSelect(1:3,3),sameline_nospace=.true.)
          if (nmsel > 3) then
             call iw_text("4",rgb=ColorMeasureSelect(1:3,4),sameline_nospace=.true.)
          else
             call iw_text("*",sameline_nospace=.true.)
          end if
          call iw_text(")=" // string(ang,'f',decimal=2) // "°",sameline_nospace=.true.)
       end if

       ! dihedral 1-2-3-4
       idx4 = msel(1:4,1)
       ang = sys(w%view_selected)%c%dihedral(&
          sys(w%view_selected)%c%atcel(idx4(1))%x + idx4(2:4),&
          sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4),&
          sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4),&
          sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4)) * 180d0 / pi
       call iw_text(", φ(",sameline_nospace=.true.)
       call iw_text("1",rgb=ColorMeasureSelect(1:3,1),sameline_nospace=.true.)
       call iw_text("2",rgb=ColorMeasureSelect(1:3,2),sameline_nospace=.true.)
       call iw_text("3",rgb=ColorMeasureSelect(1:3,3),sameline_nospace=.true.)
       if (nmsel > 3) then
          call iw_text("4",rgb=ColorMeasureSelect(1:3,4),sameline_nospace=.true.)
       else
          call iw_text("*",sameline_nospace=.true.)
       end if
       call iw_text(")=" // string(ang,'f',decimal=2) // "°",sameline_nospace=.true.)
    end if

999 continue ! exit here

    ! finish tooltip
    call igPopTextWrapPos()
    call igEndTooltip()

  end subroutine draw_selection_tooltip

  !> Mouse position to texture position (screen coordinates)
  module subroutine mousepos_to_texpos(w,pos)
    class(window), intent(inout), target :: w
    type(ImVec2), intent(inout) :: pos

    real(c_float) :: dx, dy, xratio, yratio

    if (abs(pos%x) > 1e20 .or. abs(pos%y) > 1e20) return

    dx = max(w%v_rmax%x - w%v_rmin%x,1._c_float)
    dy = max(w%v_rmax%y - w%v_rmin%y,1._c_float)
    xratio = 2._c_float * dx / max(dx,dy)
    yratio = 2._c_float * dy / max(dx,dy)

    pos%x = ((pos%x - w%v_rmin%x) / dx - 0.5_c_float) * xratio
    pos%y = (0.5_c_float - (pos%y - w%v_rmin%y) / dy) * yratio

    pos%x = (0.5_c_float * pos%x + 0.5_c_float) * w%FBOside
    ! the render texture is now presented right-side up (see draw_view), so the
    ! texture row increases upward like the framebuffer (origin bottom-left):
    ! screen-top maps to the top framebuffer row, matching glReadPixels in getpixel.
    pos%y = (0.5_c_float + 0.5_c_float * pos%y) * w%FBOside

  end subroutine mousepos_to_texpos

  !> Texture position (screen coordinates) to mouse position
  module subroutine texpos_to_mousepos(w,pos)
    class(window), intent(inout), target :: w
    type(ImVec2), intent(inout) :: pos

    real(c_float) :: x, y, xratio1, yratio1

    pos%x = (pos%x / w%FBOside) * 2._c_float - 1._c_float
    ! inverse of mousepos_to_texpos: the texture is presented right-side up, so
    ! the vertical mapping is no longer inverted
    pos%y = (pos%y / w%FBOside) * 2._c_float - 1._c_float

    x = max(w%v_rmax%x - w%v_rmin%x,1._c_float)
    y = max(w%v_rmax%y - w%v_rmin%y,1._c_float)
    xratio1 = 0.5_c_float * max(x,y) / x
    yratio1 = 0.5_c_float * max(x,y) / y

    pos%x = w%v_rmin%x + x * (0.5_c_float + xratio1 * pos%x)
    pos%y = w%v_rmin%y + y * (0.5_c_float + yratio1 * pos%y)

  end subroutine texpos_to_mousepos

  !> Get the view depth from the texture position
  module subroutine getpixel(w,pos,depth,rgba)
    use interfaces_opengl3
    class(window), intent(inout), target :: w
    type(ImVec2), intent(inout) :: pos
    real(c_float), intent(out), optional :: depth
    real(c_float), intent(out), optional :: rgba(4)

    real(c_float), target :: depth_, rgba_(4)

    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBOpick)
    if (present(depth)) then
       call glReadPixels(int(pos%x), int(pos%y), 1_c_int, 1_c_int, GL_DEPTH_COMPONENT, GL_FLOAT, c_loc(depth_))
       depth = depth_
    end if
    if (present(rgba)) then
       call glReadPixels(int(pos%x), int(pos%y), 1_c_int, 1_c_int, GL_RGBA, GL_FLOAT, c_loc(rgba_))
       rgba = rgba_
    end if
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

  end subroutine getpixel

  !> Transform from view coordinates to texture position (x,y)
  !> plus depth (z).
  module subroutine view_to_texpos(w,pos)
    use utils, only: project, unproject
    use scenes, only: scene
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    if (.not.associated(w%sc)) return
    if (w%sc%isinit < 2) return

    call project(pos,eye4,w%sc%projection,w%FBOside)

  end subroutine view_to_texpos

  !> Transform texture position (x,y) plus depth (z) to view
  !> coordinates.
  module subroutine texpos_to_view(w,pos)
    use utils, only: unproject
    use scenes, only: scene
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    if (.not.associated(w%sc)) return
    if (w%sc%isinit < 2) return

    call unproject(pos,eye4,w%sc%projection,w%FBOside)

  end subroutine texpos_to_view

  !> Transform from world coordinates to texture position (x,y)
  !> plus depth (z).
  module subroutine world_to_texpos(w,pos)
    use utils, only: project
    use scenes, only: scene
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    if (.not.associated(w%sc)) return
    if (w%sc%isinit < 2) return

    call project(pos,w%sc%view,w%sc%projection,w%FBOside)

  end subroutine world_to_texpos

  !> Transform texture position (x,y) plus depth (z) to world
  !> coordinates.
  module subroutine texpos_to_world(w,pos)
    use utils, only: unproject
    use scenes, only: scene
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    if (.not.associated(w%sc)) return
    if (w%sc%isinit < 2) return

    call unproject(pos,w%sc%view,w%sc%projection,w%FBOside)

  end subroutine texpos_to_world

  !> Export the current view to an image file with currently selected
  !> window options. The file name is file and the file format (PNG,
  !> BMP, TGA, JPE) is fformat. nsample = number of samples for
  !> anti-aliasing jpgquality = jpg quality. exportview = export the
  !> view or the whole texture.  npixel = number of pixels in the
  !> export buffer. transparentbg = transparent background. Sets the
  !> w%errmsg flag, which becomes non-zero in case of error.
  module subroutine export_to_image(w,file,fformat,nsample,npixel,transparentbg,&
     exportview,jpgquality,errmsg)
    use interfaces_opengl3
    use interfaces_stb
    use systems, only: sys_init, ok_system
    use tools_io, only: string
    class(window), intent(inout), target :: w
    character(kind=c_char,len=*), intent(in) :: file
    character(kind=c_char,len=*), intent(in) :: fformat
    integer(c_int), intent(in) :: nsample, jpgquality, npixel
    logical, intent(in) :: exportview, transparentbg
    character(kind=c_char,len=:), allocatable, intent(inout) :: errmsg

    integer(c_int), target :: msFBO, endFBO ! framebuffer
    integer(c_int), target :: msFBOdepth, endFBOdepth ! framebuffer, depth buffer
    integer(c_int), target :: msFBOtex, endFBOtex ! framebuffer, texture
    integer :: isys
    integer :: width, height, origin(2)
    integer(c_signed_char), allocatable, target :: data(:)
    type(ImVec2) :: x0, x1
    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: idum

    ! reset the error message
    errmsg = ""

    ! get the scene and the system
    if (.not.associated(w%sc)) then
       errmsg = "Scene is not available for exporting"
       return
    end if
    isys = w%sc%id
    if (.not.ok_system(isys,sys_init) .or. isys /= w%view_selected) then
       errmsg = "System is not initialized"
       return
    end if

    ! generate textures and buffers
    call glGenTextures(1, c_loc(msFBOtex))
    call glGenRenderbuffers(1, c_loc(msFBOdepth))
    call glGenFramebuffers(1, c_loc(msFBO))
    call glGenTextures(1, c_loc(endFBOtex))
    call glGenRenderbuffers(1, c_loc(endFBOdepth))
    call glGenFramebuffers(1, c_loc(endFBO))

    ! textures
    call glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, msFBOtex)
    call glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, nsample, GL_RGBA, npixel, npixel,&
       int(GL_TRUE,c_signed_char))
    call glBindTexture(GL_TEXTURE_2D, 0)
    call glBindTexture(GL_TEXTURE_2D, endFBOtex)
    call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, npixel, npixel, 0, GL_RGBA, GL_UNSIGNED_BYTE, c_null_ptr)
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    call glBindTexture(GL_TEXTURE_2D, 0)

    ! render buffer
    call glBindRenderbuffer(GL_RENDERBUFFER, msFBOdepth)
    call glRenderbufferStorageMultisample(GL_RENDERBUFFER, nsample, GL_DEPTH_COMPONENT, npixel, npixel)
    call glBindRenderbuffer(GL_RENDERBUFFER, 0)
    call glBindRenderbuffer(GL_RENDERBUFFER, endFBOdepth)
    call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, npixel, npixel)
    call glBindRenderbuffer(GL_RENDERBUFFER, 0)

    ! frame buffer
    call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
    call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, msFBOtex, 0)
    call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, msFBOdepth)
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
       errmsg = "Failed creating the multi-sampled framebuffer (too large?)"
       goto 999
    end if
    call glBindFramebuffer(GL_FRAMEBUFFER, endFBO)
    call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, endFBOtex, 0)
    call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, endFBOdepth)
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
       errmsg = "Failed creating the render framebuffer (too large?)"
       goto 999
    end if
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

    ! render the scene to the multisampled framebuffer
    call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
    call glViewport(0_c_int,0_c_int,npixel,npixel)
    if (transparentbg) then
       call glClearColor(w%sc%bgcolor(1),w%sc%bgcolor(2),w%sc%bgcolor(3),0._c_float)
    else
       call glClearColor(w%sc%bgcolor(1),w%sc%bgcolor(2),w%sc%bgcolor(3),1._c_float)
    end if
    call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
    call w%sc%render()
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

    ! blit the multisampled buffer to the normal colorbuffer
    call glBindFramebuffer(GL_READ_FRAMEBUFFER, msFBO)
    call glBindFramebuffer(GL_DRAW_FRAMEBUFFER, endFBO)
    call glBlitFramebuffer(0, 0, npixel, npixel, 0, 0, npixel, npixel, GL_COLOR_BUFFER_BIT, GL_LINEAR)
    call glBindFramebuffer(GL_READ_FRAMEBUFFER, 0)
    call glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0)

    ! Read from the regular framebuffer into the data array
    call glBindFramebuffer(GL_FRAMEBUFFER, endFBO)
    if (.not.exportview) then
       ! whole texture
       width = npixel
       height = npixel
       allocate(data(4 * width * height))
       data = 0_c_signed_char
       call glReadPixels(0, 0, npixel, npixel, GL_RGBA, GL_UNSIGNED_BYTE, c_loc(data))
    else
       ! viewport only
       x0%x = w%v_rmin%x
       x0%y = w%v_rmin%y
       x1%x = w%v_rmax%x
       x1%y = w%v_rmax%y
       call w%mousepos_to_texpos(x0)
       call w%mousepos_to_texpos(x1)
       ! the texture is presented right-side up, so v_rmin (screen top) maps to
       ! the larger texture row and v_rmax to the smaller; use abs/min so the
       ! crop rectangle (glReadPixels origin is the lower-left) is well defined
       width = min(nint(abs(x1%x - x0%x) / real(w%FBOside,8) * npixel),npixel)
       height = min(nint(abs(x1%y - x0%y) / real(w%FBOside,8) * npixel),npixel)
       origin(1) = max(nint(min(x0%x,x1%x) / real(w%FBOside,8) * npixel),0)
       origin(2) = max(nint(min(x0%y,x1%y) / real(w%FBOside,8) * npixel),0)
       allocate(data(4 * width * height))
       data = 0_c_signed_char
       call glReadPixels(origin(1), origin(2), width, height, GL_RGBA, GL_UNSIGNED_BYTE, c_loc(data))
    end if
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
       errmsg = "Error rendering export image"

    ! write the file. glReadPixels returns rows bottom-to-top while stb writes
    ! top-to-bottom, so flip vertically on write to match the on-screen view.
    call stbi_flip_vertically_on_write(1_c_int)
    str = trim(file) // c_null_char
    if (fformat(1:3) == "PNG") then
       idum = stbi_write_png(c_loc(str), width, height, 4, c_loc(data), 4*width)
    elseif (fformat(1:3) == "BMP") then
       idum = stbi_write_bmp(c_loc(str), width, height, 4, c_loc(data))
    elseif (fformat(1:3) == "TGA") then
       idum = stbi_write_tga(c_loc(str), width, height, 4, c_loc(data))
    elseif (fformat(1:3) == "JPE") then
       idum = stbi_write_jpg(c_loc(str), width, height, 4, c_loc(data), jpgquality)
    else
       idum = 0
       errmsg = "Unknown file format: "  // string(fformat)
    end if
    call stbi_flip_vertically_on_write(0_c_int)
    if (idum == 0) &
       errmsg = "Error exporting image to file: "  // string(file)

999 continue

    ! delete the buffers
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)
    call glDeleteTextures(1, c_loc(msFBOtex))
    call glDeleteRenderbuffers(1, c_loc(msFBOdepth))
    call glDeleteFramebuffers(1, c_loc(msFBO))
    call glDeleteTextures(1, c_loc(endFBOtex))
    call glDeleteRenderbuffers(1, c_loc(endFBOdepth))
    call glDeleteFramebuffers(1, c_loc(endFBO))

  end subroutine export_to_image

  !> Export the current system to a PNG file using the same file name
  !> but with png extension, using all defaults
  module subroutine export_to_png_simple(w)
    use tools_io, only: equal, lower
    use systems, only: sys_loaded_not_init, sysc, ok_system
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable :: file, errmsg
    character(len=:), allocatable :: wext, wroot
    integer(c_int) :: npixel
    integer :: isys, idx

    character(kind=c_char,len=3), parameter :: fformat = "PNG"
    integer(c_int), parameter :: nsample = 16
    integer(c_int), parameter :: jpgquality = 90
    logical, parameter :: exportview = .true.
    logical, parameter :: transparentbg = .true.

    ! initialize
    npixel = w%FBOside

    ! get the root of the file name and append the png extension
    isys = w%view_selected
    if (.not.ok_system(isys,sys_loaded_not_init)) return
    file = sysc(isys)%seed%file
    wext = lower(file(index(file,'.',.true.)+1:))
    wroot = file(:index(file,'.',.true.)-1)
    if (equal(wext,'in')) then
       idx = index(wroot,'.',.true.)
       if (idx > 0) then
          wext = lower(wroot(idx+1:))
          if (equal(wext,'scf')) &
             wroot = file(:index(file,'.',.true.)-1)
       end if
    end if
    file = wroot // ".png"

    ! export
    call w%export_to_image(file,fformat,nsample,npixel,transparentbg,exportview,&
       jpgquality,errmsg)

  end subroutine export_to_png_simple

end submodule view
