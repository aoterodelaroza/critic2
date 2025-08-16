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

contains

  !> Draw the view.
  module subroutine draw_view(w)
    use interfaces_opengl3
    use interfaces_cimgui
    use keybindings, only: is_bind_event, BIND_VIEW_INC_NCELL, BIND_VIEW_DEC_NCELL,&
       BIND_VIEW_ALIGN_A_AXIS, BIND_VIEW_ALIGN_B_AXIS, BIND_VIEW_ALIGN_C_AXIS,&
       BIND_VIEW_ALIGN_X_AXIS, BIND_VIEW_ALIGN_Y_AXIS, BIND_VIEW_ALIGN_Z_AXIS,&
       BIND_VIEW_TOGGLE_ATOMS, BIND_VIEW_TOGGLE_BONDS, BIND_VIEW_CYCLE_LABELS,&
       BIND_VIEW_TOGGLE_CELL,&
       BIND_NAV_ROTATE, BIND_NAV_ROTATE_PERP, BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET,&
       BIND_NAV_MEASURE, bindnames, get_bind_keyname,&
       BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use representations, only: reptype_atoms, reptype_unitcell,&
       repflavor_atoms_ballandstick, repflavor_atoms_vdwcontacts, repflavor_atoms_hbonds,&
       repflavor_atoms_sticks, repflavor_atoms_licorice, repflavor_unitcell_basic
    use scenes, only: style_phong, style_simple
    use utils, only: iw_calcheight, iw_calcwidth, iw_clamp_color3, iw_combo_simple,&
       iw_setposx_fromend, iw_checkbox, iw_coloredit, iw_menuitem
    use crystalmod, only: iperiod_vacthr
    use global, only: dunit0, iunit_ang
    use systems, only: sysc, sys, sys_init, nsys, ok_system, are_threads_running
    use gui_main, only: g, fontsize, lockbehavior, tree_select_updates_view
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    integer :: i, j, k, nrep, ipad, is, icel, ineq, iaux
    type(ImVec2) :: szavail, sz0, sz1, szero, pos
    type(ImVec4) :: tintcol, bgcol
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    character(len=:), allocatable, target :: msg
    logical(c_bool) :: is_selected
    logical :: hover, chbuild, chrender, goodsys, ldum, ok, ismol, isatom, isbond
    logical :: isuc, islabelsl, hover_and_moved, enabled
    integer :: islabels
    logical(c_bool) :: ch
    integer(c_int) :: flags, nc(3), ires, viewtype, idum
    real(c_float) :: scal, width, sqw, depth, rgba(4)
    real*8 :: x0(3)
    logical :: changedisplay(4) ! 1=atoms, 2=bonds, 3=labels, 4=cell

    logical, save :: ttshown = .false. ! tooltip flag

    ! coordinate this with objects (representation) menu in scenes module
    integer(c_int), parameter :: ic_closebutton = 0
    integer(c_int), parameter :: ic_viewbutton = 1
    integer(c_int), parameter :: ic_name = 2
    integer(c_int), parameter :: ic_type = 3
    integer(c_int), parameter :: ic_editbutton = 4

    ! initialize
    szero%x = 0
    szero%y = 0
    chrender = .false.
    chbuild = .false.
    if (w%firstpass) then
       w%mousepos_lastframe%x = 0._c_float
       w%mousepos_lastframe%y = 0._c_float
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
    if (associated(w%sc)) then
       do i = 1, w%sc%nrep
          if (w%sc%rep(i)%isinit) then
             if (w%sc%rep(i)%type == reptype_atoms) then
                isatom = isatom .or. w%sc%rep(i)%atoms_display
                isbond = isbond .or. w%sc%rep(i)%bonds_display
                if (w%sc%rep(i)%labels_display) &
                   islabels = w%sc%rep(i)%label_type
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
       if (any(changedisplay)) then
          do i = 1, w%sc%nrep
             if (w%sc%rep(i)%isinit) then
                if (w%sc%rep(i)%type == reptype_atoms) then
                   if (changedisplay(1)) w%sc%rep(i)%atoms_display = isatom
                   if (changedisplay(2)) w%sc%rep(i)%bonds_display = isbond
                   if (changedisplay(3)) then
                      w%sc%rep(i)%labels_display = islabelsl
                      if (islabelsl) then
                         w%sc%rep(i)%label_type = islabels
                         call w%sc%rep(i)%label_style%reset(w%sc%rep(i))
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
                      w%sc%rep(i)%atoms_display = isatom
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
                      w%sc%rep(i)%bonds_display = isbond
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
                      w%sc%rep(i)%labels_display = islabelsl
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

             ! number of cells in each direction
             ! calculate widths
             ipad = ceiling(log10(max(maxval(w%sc%nc),1) + 0.1))
             sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
             call igPushItemWidth(sqw)

             nc = w%sc%nc
             call igAlignTextToFramePadding()
             call iw_text("a:")
             call igSameLine(0._c_float,0._c_float)
             if (iw_button("-##aaxis")) nc(1) = max(nc(1)-1,1)
             call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
             str2 = "##aaxis" // c_null_char
             call igPushItemWidth(iw_calcwidth(ipad,1))
             ldum = igInputInt(c_loc(str2),nc(1),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
             call igPopItemWidth()
             call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
             if (iw_button("+##aaxis")) nc(1) = nc(1)+1

             call igSameLine(0._c_float,-1._c_float)
             call iw_text("b:")
             call igSameLine(0._c_float,0._c_float)
             if (iw_button("-##baxis")) nc(2) = max(nc(2)-1,1)
             call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
             str2 = "##baxis" // c_null_char
             call igPushItemWidth(iw_calcwidth(ipad,1))
             ldum = igInputInt(c_loc(str2),nc(2),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
             call igPopItemWidth()
             call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
             if (iw_button("+##baxis")) nc(2) = nc(2)+1

             call igSameLine(0._c_float,-1._c_float)
             call iw_text("c:")
             call igSameLine(0._c_float,0._c_float)
             if (iw_button("-##caxis")) nc(3) = max(nc(3)-1,1)
             call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
             str2 = "##caxis" // c_null_char
             call igPushItemWidth(iw_calcwidth(ipad,1))
             ldum = igInputInt(c_loc(str2),nc(3),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
             call igPopItemWidth()
             call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
             if (iw_button("+##caxis")) nc(3) = nc(3)+1

             nc = max(nc,1)
             if (any(nc /= w%sc%nc)) then
                w%sc%nc = nc
                chbuild = .true.
             end if
             call igPopItemWidth()
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
          call igSameLine(0._c_float,-1._c_float)
          str2 = "Reset Distance##resetdistance" // c_null_char
          str3 = "%.2f" // c_null_char
          call igPushItemWidth(iw_calcwidth(5,1))
          ch = igDragFloat(c_loc(str2),w%sc%camresetdist,&
             0.01_c_float,0.1_c_float,8.0_c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
          call igPopItemWidth()
          call iw_tooltip("Ratio controlling distance from object when resetting camera",ttshown)

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

          ! scene style
          call iw_text("Appearance",highlight=.true.)
          call iw_combo_simple("Style##scenestyle","Simple"//c_null_char//"Realistic"&
             //c_null_char//c_null_char,w%sc%style,changed=ch)
          if (ch) then
             call w%sc%set_style_defaults()
             w%sc%forcebuildlists = .true.
          end if

          if (w%sc%style == style_phong) then
             !! phong-specific options !!
             call igPushItemWidth(iw_calcwidth(15,3))
             str2 = "Light Position" // c_null_char
             str3 = "%.1f" // c_null_char
             chrender = chrender .or. igDragFloat3(c_loc(str2),w%sc%lightpos,&
                0.5_c_float,-FLT_MAX,FLT_MAX,c_loc(str3),ImGuiSliderFlags_None)
             call iw_tooltip("Change the position of the light",ttshown)
             call igPopItemWidth()

             call igPushItemWidth(iw_calcwidth(5,1))
             str2 = "Ambient " // c_null_char
             str3 = "%.3f" // c_null_char
             chrender = chrender .or. igDragFloat(c_loc(str2),w%sc%ambient,&
                0.002_c_float,0._c_float,1._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the ambient light intensity",ttshown)
             call igSameLine(0._c_float,-1._c_float)
             str2 = "Diffuse" // c_null_char
             str3 = "%.3f" // c_null_char
             chrender = chrender .or. igDragFloat(c_loc(str2),w%sc%diffuse,&
                0.002_c_float,0._c_float,1._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the diffuse light intensity",ttshown)
             str2 = "Specular" // c_null_char
             str3 = "%.3f" // c_null_char
             chrender = chrender .or. igDragFloat(c_loc(str2),w%sc%specular,&
                0.002_c_float,0._c_float,1._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the specular light intensity",ttshown)
             call igSameLine(0._c_float,-1._c_float)
             str2 = "Shininess" // c_null_char
             str3 = "%.0f" // c_null_char
             chrender = chrender .or. igDragInt(c_loc(str2),w%sc%shininess,&
                1._c_float,0_c_int,256_c_int,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the shininess of the light",ttshown)
             call igPopItemWidth()

             chrender = chrender .or. iw_coloredit("Light",rgb=w%sc%lightcolor)
             call iw_tooltip("Change the color of the light",ttshown)
             call igSameLine(0._c_float,-1._c_float)
          elseif (w%sc%style == style_simple) then
             call igAlignTextToFramePadding()
             call iw_text("Atom Border: ")
             call igSameLine(0._c_float,-1._c_float)

             chrender = chrender .or. iw_coloredit("Color",rgb=w%sc%bordercolor,sameline=.true.)
             call iw_tooltip("Change the color of the atom borders",ttshown)
          end if

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
                            sysc(i)%sc%rep(j)%atoms_display = isatom
                            sysc(i)%sc%rep(j)%bonds_display = isbond
                            sysc(i)%sc%rep(j)%labels_display = islabelsl
                            if (islabelsl) then
                               if (sys(i)%c%ismolecule.and.islabels == 8) then
                                  sysc(i)%sc%rep(j)%label_type = 0
                               else
                                  sysc(i)%sc%rep(j)%label_type = islabels
                               end if
                               call sysc(i)%sc%rep(j)%label_style%reset(sysc(i)%sc%rep(j))
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
                   sysc(i)%sc%lightpos = w%sc%lightpos
                   sysc(i)%sc%ambient = w%sc%ambient
                   sysc(i)%sc%diffuse = w%sc%diffuse
                   sysc(i)%sc%specular = w%sc%specular
                   sysc(i)%sc%shininess = w%sc%shininess
                   sysc(i)%sc%bordercolor = w%sc%bordercolor
                   sysc(i)%sc%lightcolor = w%sc%lightcolor
                   sysc(i)%sc%bgcolor = w%sc%bgcolor
                   sysc(i)%sc%camresetdist = w%sc%camresetdist
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

             if (.not.sys(w%view_selected)%c%ismolecule) then
                call igSeparator()
                if (iw_menuitem("Unit Cell",shortcut_text="Cell")) then
                   call w%sc%add_representation(reptype_unitcell,repflavor_unitcell_basic)
                   chbuild = .true.
                end if
                call iw_tooltip("Display the unit cell",ttshown)
             end if

             call igEndPopup()
          end if
          call iw_tooltip("Add a new object to the view",ttshown)

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

          call igEndPopup()
       end if
    end if
    call iw_tooltip("Add, remove, and modify objects",ttshown)

    ! update the draw lists and render
    if (associated(w%sc)) then
       if (chbuild .or. w%sc%timelastbuild < sysc(w%view_selected)%timelastchange_buildlists) &
          w%sc%forcebuildlists = .true.
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

       if (iw_menuitem("Recalculate Bonds...",enabled=enabled.and..not.are_threads_running())) &
          idum = stack_create_window(wintype_rebond,.true.,isys=w%view_selected,orraise=-1)
       call iw_tooltip("Recalculate the bonds in the current system",ttshown)

       if (iw_menuitem("Vibrations...",enabled=enabled)) &
          iaux = stack_create_window(wintype_vibrations,.true.,idparent=w%id,orraise=-1)
       call iw_tooltip("Display an animation showing the atomic vibrations for this system",ttshown)
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

    ! ! resize the render texture if not large enough
    !!! integer(c_int) :: amax
    ! amax = max(ceiling(max(szavail%x,szavail%y)),1)
    ! if (amax > w%FBOside) then
    !    amax = max(ceiling(1.2 * ceiling(max(szavail%x,szavail%y))),1)
    !    call w%delete_texture_view()
    !    call w%create_texture_view(amax)
    !    w%forcerender = .true.
    ! end if

    ! draw the texture, largest region with the same shape as the available region
    ! that fits into the texture square
    scal = real(w%FBOside,c_float) / max(max(szavail%x,szavail%y),1._c_float)
    sz0%x = 0.5 * (real(w%FBOside,c_float) - szavail%x * scal) / real(w%FBOside,c_float)
    sz0%y = 0.5 * (real(w%FBOside,c_float) - szavail%y * scal) / real(w%FBOside,c_float)
    sz1%x = 1._c_float - sz0%x
    sz1%y = 1._c_float - sz0%y

    !!! The camratio is used for resetting the camera. The problem with this
    !!! is that in the first render the window dimensions are still changing
    !!! and the camratio of the first view rendered is different from the rest.
    !!! Therefore, set to a constant of 1.5 in scene_init for the time being.
    ! if (szavail%x > szavail%y) then
    !    ratio = szavail%x / max(szavail%y,1._c_float)
    ! else
    !    ratio = szavail%y / max(szavail%x,1._c_float)
    ! end if
    ! if (associated(w%sc)) then
    !    w%sc%camratio = min(ratio,2.5_c_float)
    ! end if

    ! render the image to the texture, if requested
    if (w%forcerender) then
       ! render to the draw framebuffer
       call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
       call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
       if (associated(w%sc)) then
          call glClearColor(w%sc%bgcolor(1),w%sc%bgcolor(2),&
             w%sc%bgcolor(3),1._c_float)
       else
          call glClearColor(0._c_float,0._c_float,0._c_float,0._c_float)
       end if
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
       if (associated(w%sc)) call w%sc%render()
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! render to the pick frame buffer
       call glBindFramebuffer(GL_FRAMEBUFFER, w%FBOpick)
       call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
       call glClearColor(0._c_float,0._c_float,0._c_float,0._c_float)
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
       if (associated(w%sc)) call w%sc%renderpick()
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       w%forcerender = .false.
    end if

    ! border and tint for the image, draw the image, update the rectangle
    tintcol%x = 1._c_float
    tintcol%y = 1._c_float
    tintcol%z = 1._c_float
    tintcol%w = 1._c_float
    bgcol%x = 0._c_float
    bgcol%y = 0._c_float
    bgcol%z = 0._c_float
    bgcol%w = 1._c_float
    call igPushStyleColor_Vec4(ImGuiCol_Button,bgcol)
    call igPushStyleColor_Vec4(ImGuiCol_ButtonActive,bgcol)
    call igPushStyleColor_Vec4(ImGuiCol_ButtonHovered,bgcol)
    str1 = "##imagebutton" // c_null_char
    ldum = igImageButtonEx(igGetID_Str(c_loc(str1)),w%FBOtex, szavail, sz0, sz1, szero, bgcol, tintcol)
    call igPopStyleColor(3)

    ! get hover, hover_and_moved, image rectangle coordinates, and atom idx
    hover = goodsys .and. w%ilock == ilock_no
    if (hover) hover = igIsItemHovered(ImGuiHoveredFlags_None)
    hover_and_moved = .false.
    if (hover) then
       call igGetMousePos(pos)
       hover_and_moved = (abs(w%mousepos_lastframe%x-pos%x) > 1e-4.or.abs(w%mousepos_lastframe%y-pos%y) > 1e-4)
       hover_and_moved = hover_and_moved.or.(w%view_mousebehavior == MB_Navigation.and.is_bind_event(BIND_NAV_MEASURE))
       w%mousepos_lastframe = pos
    end if

    ! get the ID of the atom under mouse
    call igGetItemRectMin(w%v_rmin)
    call igGetItemRectMax(w%v_rmax)
    if (hover_and_moved) then
       w%mousepos_idx = 0
       call w%mousepos_to_texpos(pos)
       call w%getpixel(pos,depth,rgba)

       ! transform to atom cell ID and lattice vector
       w%mousepos_idx(1:4) = transfer(rgba,w%mousepos_idx(1:4))
       w%mousepos_idx(5) = w%mousepos_idx(1)
       if (associated(w%sc) .and. w%mousepos_idx(1) > 0 .and. w%mousepos_idx(1) <= w%sc%obj%nsph) then
          w%mousepos_idx(1:4) = w%sc%obj%sph(w%mousepos_idx(1))%idx
       else
          w%mousepos_idx = 0
       end if
    elseif (.not.hover) then
       w%mousepos_idx = 0
    end if

    ! mode selection
    viewtype = 0
    call iw_combo_simple("##viewmode","Navigate" // c_null_char,viewtype)
    msg = trim(get_bind_keyname(BIND_NAV_ROTATE)) // ": " // trim(bindnames(BIND_NAV_ROTATE)) // newline
    msg = msg // trim(get_bind_keyname(BIND_NAV_ROTATE_PERP)) // ": " // trim(bindnames(BIND_NAV_ROTATE_PERP)) // newline
    msg = msg // trim(get_bind_keyname(BIND_NAV_TRANSLATE)) // ": " // trim(bindnames(BIND_NAV_TRANSLATE)) // newline
    msg = msg // trim(get_bind_keyname(BIND_NAV_ZOOM)) // ": " // trim(bindnames(BIND_NAV_ZOOM)) // newline
    msg = msg // trim(get_bind_keyname(BIND_NAV_RESET)) // ": " // trim(bindnames(BIND_NAV_RESET)) // newline
    msg = msg // trim(get_bind_keyname(BIND_NAV_MEASURE)) // ": " // trim(bindnames(BIND_NAV_MEASURE)) // newline
    call iw_tooltip(msg)

    ! atom hover message
    if (hover .and. w%mousepos_idx(1) > 0) then
       call igSameLine(0._c_float,-1._c_float)
       icel = w%mousepos_idx(1)
       is = sys(w%view_selected)%c%atcel(icel)%is
       ineq = sys(w%view_selected)%c%atcel(icel)%idx
       ismol = sys(w%view_selected)%c%ismolecule

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
    call w%process_events_view(hover,w%mousepos_idx)

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

    if (.not.w%ismain) then
       ! the close button
       if (iw_button("Close",danger=.true.)) w%isopen = .false.

       ! exit if focused and received the close keybinding
       if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) .or.&
          is_bind_event(BIND_CLOSE_ALL_DIALOGS)) then
          w%isopen = .false.
       end if
    end if

  end subroutine draw_view

  !> Create the texture for the view window, with atex x atex pixels.
  module subroutine create_texture_view(w,atex)
    use interfaces_opengl3
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
    call glClearColor(0._c_float,0._c_float,0._c_float,0._c_float)
    call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)

    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBOpick)
    call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
    call glClearColor(0._c_float,0._c_float,0._c_float,0._c_float)
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
    w%mousepos_lastframe%x = 0._c_float
    w%mousepos_lastframe%y = 0._c_float
    w%mousepos_idx = 0

    ! set the time
    w%timelast_view_assign = glfwGetTime()

  end subroutine select_view

  !> Process the mouse events in the view window
  module subroutine process_events_view(w,hover,idx)
    use interfaces_cimgui
    use scenes, only: scene
    use utils, only: translate, rotate, mult, invmult
    use tools_math, only: cross_cfloat, matinv_cfloat
    use keybindings, only: is_bind_event, is_bind_mousescroll, BIND_NAV_ROTATE,&
       BIND_NAV_ROTATE_PERP,&
       BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET, BIND_NAV_MEASURE,&
       BIND_CLOSE_FOCUSED_DIALOG
    use systems, only: nsys
    use gui_main, only: io
    class(window), intent(inout), target :: w
    logical, intent(in) :: hover
    integer(c_int), intent(in) :: idx(5)

    type(ImVec2) :: texpos, mousepos
    real(c_float) :: ratio, pos3(3), vnew(3), vold(3), axis(3)
    real(c_float) :: mpos2(2), ang, xc(3), dist

    real(c_float), parameter :: mousesens_zoom0 = 0.15_c_float
    real(c_float), parameter :: mousesens_rot0 = 3._c_float

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
    end if

    ! only process if there is an associated system is viewed and scene is initialized
    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    if (.not.associated(w%sc)) return
    if (w%sc%isinit < 2) return

    ! process mode-specific events
    if (w%view_mousebehavior == MB_Navigation) then
       call igGetMousePos(mousepos)
       texpos = mousepos

       ! transform to the texture pos
       call w%mousepos_to_texpos(texpos)

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
          w%mpos0_r = (/texpos%x,texpos%y,1._c_float/)

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
          call w%sc%select_atom(idx)
          w%forcerender = .true.
       end if
       if (hover .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) then
          call w%sc%select_atom((/0,0,0,0,0/))
          w%forcerender = .true.
       end if
    end if

  end subroutine process_events_view

  module subroutine draw_selection_tooltip(w,idx)
    use interfaces_cimgui
    use utils, only: iw_text
    use systems, only: sys
    use gui_main, only: fontsize, ColorMeasureSelect
    use tools_io, only: string
    use tools_math, only: cross
    use param, only: bohrtoa, pi
    class(window), intent(inout), target :: w
    integer(c_int), intent(in) :: idx(5)

    integer :: nmsel
    integer :: msel(5,4)
    integer :: idx1(4), idx2(4), idx3(4), idx4(4)
    real*8 :: x0(3), x1(3), x2(3), d, d1, d2, ang, n0(3), n1(3)

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
       x0 = sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4) -&
          (sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       x1 = sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4) -&
          (sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       x0 = sys(w%view_selected)%c%x2c(x0)
       x1 = sys(w%view_selected)%c%x2c(x1)
       d1 = norm2(x0)
       d2 = norm2(x1)
       if (d1 > 1d-14 .and. d2 > 1d-14) then
          ang = acos(dot_product(x0,x1) / d1 / d2) * 180d0 / pi
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
       x0 = sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4) -&
          (sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       x1 = sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4) -&
          (sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       x0 = sys(w%view_selected)%c%x2c(x0)
       x1 = sys(w%view_selected)%c%x2c(x1)
       d1 = norm2(x0)
       d2 = norm2(x1)
       if (d1 > 1d-14 .and. d2 > 1d-14) then
          ang = acos(dot_product(x0,x1) / norm2(x0) / norm2(x1)) * 180d0 / pi
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
       x0 = sys(w%view_selected)%c%atcel(idx4(1))%x + idx4(2:4) -&
          (sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4))
       x1 = sys(w%view_selected)%c%atcel(idx3(1))%x + idx3(2:4) -&
          (sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4))
       x2 = sys(w%view_selected)%c%atcel(idx1(1))%x + idx1(2:4) -&
          (sys(w%view_selected)%c%atcel(idx2(1))%x + idx2(2:4))
       x0 = sys(w%view_selected)%c%x2c(x0)
       x1 = sys(w%view_selected)%c%x2c(x1)
       x2 = sys(w%view_selected)%c%x2c(x2)
       n0 = cross(x0,x1)
       n1 = cross(x1,x2)

       ang = -atan2(norm2(x1) * dot_product(x0,n1), dot_product(n0,n1)) * 180d0/pi
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
    pos%y = (0.5_c_float - 0.5_c_float * pos%y) * w%FBOside

  end subroutine mousepos_to_texpos

  !> Texture position (screen coordinates) to mouse position
  module subroutine texpos_to_mousepos(w,pos)
    class(window), intent(inout), target :: w
    type(ImVec2), intent(inout) :: pos

    real(c_float) :: x, y, xratio1, yratio1

    pos%x = (pos%x / w%FBOside) * 2._c_float - 1._c_float
    pos%y = 1._c_float - (pos%y / w%FBOside) * 2._c_float

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
       width = min(nint((x1%x - x0%x) / real(w%FBOside,8) * npixel),npixel)
       height = min(nint((x1%y - x0%y) / real(w%FBOside,8) * npixel),npixel)
       origin(1) = max(nint(x0%x / real(w%FBOside,8) * npixel),0)
       origin(2) = max(nint(x0%y / real(w%FBOside,8) * npixel),0)
       allocate(data(4 * width * height))
       data = 0_c_signed_char
       call glReadPixels(origin(1), origin(2), width, height, GL_RGBA, GL_UNSIGNED_BYTE, c_loc(data))
    end if
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
       errmsg = "Error rendering export image"

    ! write the file
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
