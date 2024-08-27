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

  !xx! private procedures
  ! function atom_selection_widget() result(changed)

contains

  !xx! view

  !> Draw the view.
  module subroutine draw_view(w)
    use interfaces_opengl3
    use interfaces_cimgui
    use keybindings, only: is_bind_event, BIND_VIEW_INC_NCELL, BIND_VIEW_DEC_NCELL,&
       BIND_VIEW_ALIGN_A_AXIS, BIND_VIEW_ALIGN_B_AXIS, BIND_VIEW_ALIGN_C_AXIS,&
       BIND_VIEW_ALIGN_X_AXIS, BIND_VIEW_ALIGN_Y_AXIS, BIND_VIEW_ALIGN_Z_AXIS,&
       BIND_VIEW_TOGGLE_ATOMS, BIND_VIEW_TOGGLE_BONDS, BIND_VIEW_TOGGLE_LABELS,&
       BIND_VIEW_TOGGLE_CELL,&
       BIND_NAV_ROTATE, BIND_NAV_ROTATE_PERP, BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET,&
       BIND_NAV_MEASURE, bindnames, get_bind_keyname,&
       BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use scenes, only: reptype_atoms, reptype_unitcell, style_phong, style_simple,&
       repflavor_atoms_basic, repflavor_atoms_vdwcontacts, repflavor_atoms_hbonds,&
       repflavor_unitcell_basic
    use utils, only: iw_calcheight, iw_calcwidth, iw_clamp_color3, iw_combo_simple,&
       iw_setposx_fromend, iw_checkbox, iw_coloredit3
    use crystalmod, only: iperiod_vacthr
    use global, only: dunit0, iunit_ang
    use gui_main, only: sysc, sys, sys_init, nsys, g, fontsize, lockbehavior
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    integer :: i, j, k, nrep, id, ipad, is, icel, ineq
    type(ImVec2) :: szavail, sz0, sz1, szero, pos
    type(ImVec4) :: tintcol, bgcol
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    character(len=:), allocatable, target :: msg
    logical(c_bool) :: is_selected
    logical :: hover, chbuild, chrender, goodsys, ldum, ok, ismol, isatom, isbond
    logical :: islabels, isuc
    logical(c_bool) :: ch, enabled
    integer(c_int) :: flags, nc(3), ires, idx(5), viewtype, idum
    real(c_float) :: scal, width, sqw, ratio, depth, rgba(4)
    real*8 :: x0(3)
    logical :: changedisplay(4) ! 1=atoms, 2=bonds, 3=labels, 4=cell

    logical, save :: ttshown = .false. ! tooltip flag

    ! coordinate this with objects (representation) menu in scenes module
    integer(c_int), parameter :: ic_closebutton = 0
    integer(c_int), parameter :: ic_viewbutton = 1
    integer(c_int), parameter :: ic_name = 2
    integer(c_int), parameter :: ic_editbutton = 3

    ! initialize
    szero%x = 0
    szero%y = 0
    chrender = .false.
    chbuild = .false.

    ! update ID for the export window
    call update_window_id(w%idexportwin)
    call update_window_id(w%idvibrationswin)

    ! whether the selected view system is a good system, and associate the scene
    goodsys = (w%view_selected >= 1 .and. w%view_selected <= nsys)
    if (goodsys) goodsys = sysc(w%view_selected)%status == sys_init
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
    islabels = .false.
    isuc = .false.
    if (associated(w%sc)) then
       do i = 1, w%sc%nrep
          if (w%sc%rep(i)%isinit) then
             if (w%sc%rep(i)%type == reptype_atoms) then
                isatom = isatom .or. w%sc%rep(i)%atoms_display
                isbond = isbond .or. w%sc%rep(i)%bonds_display
                islabels = islabels .or. w%sc%rep(i)%labels_display
             elseif (w%sc%rep(i)%type == reptype_unitcell) then
                isuc = isuc .or. w%sc%rep(i)%shown
             end if
          end if
       end do
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
       if (is_bind_event(BIND_VIEW_TOGGLE_LABELS)) then
          changedisplay(3) = .true.
          islabels = .not.islabels
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
                   if (changedisplay(3)) w%sc%rep(i)%labels_display = islabels
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

          if (iw_checkbox("Labels##labelshortcut",islabels,sameline=.true.)) then
             do i = 1, w%sc%nrep
                if (w%sc%rep(i)%isinit) then
                   if (w%sc%rep(i)%type == reptype_atoms) then
                      w%sc%rep(i)%labels_display = islabels
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display labels in all objects ("//&
             trim(get_bind_keyname(BIND_VIEW_TOGGLE_LABELS)) // ").",ttshown)

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

             chrender = chrender .or. iw_coloredit3("Light",w%sc%lightcolor)
             call iw_tooltip("Change the color of the light",ttshown)
             call igSameLine(0._c_float,-1._c_float)
          elseif (w%sc%style == style_simple) then
             call igAlignTextToFramePadding()
             call iw_text("Atom Border: ")
             call igSameLine(0._c_float,-1._c_float)

             call igPushItemWidth(iw_calcwidth(5,1))
             str2 = "Size (Å)" // c_null_char
             str3 = "%.3f" // c_null_char
             chrender = chrender .or. igDragFloat(c_loc(str2),w%sc%atomborder,&
                0.002_c_float,0._c_float,1._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the thickness of the atom borders",ttshown)
             call igPopItemWidth()

             chrender = chrender .or. iw_coloredit3("Color",w%sc%bordercolor,sameline=.true.)
             call iw_tooltip("Change the color of the atom borders",ttshown)
          end if

          ! background color
          chrender = chrender .or. iw_coloredit3("Background",w%sc%bgcolor)
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
                            sysc(i)%sc%rep(j)%labels_display = islabels
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
                   sysc(i)%sc%atomborder = w%sc%atomborder
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
          ! objectss table
          call iw_text("List of Objects",highlight=.true.)

          ! add button
          ldum = iw_button("Add",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
          if (ok) then
             str2 = "Atoms" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                id = w%sc%get_new_representation_id()
                call w%sc%rep(id)%init(w%sc,w%view_selected,id,reptype_atoms,w%sc%style,&
                   repflavor_atoms_basic)
                chbuild = .true.
             end if
             call iw_tooltip("Display atoms, bonds, and labels in the scene",ttshown)

             if (.not.sys(w%view_selected)%c%ismolecule) then
                str2 = "Unit Cell" // c_null_char
                if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                   id = w%sc%get_new_representation_id()
                   call w%sc%rep(id)%init(w%sc,w%view_selected,id,reptype_unitcell,w%sc%style,&
                      repflavor_unitcell_basic)
                   chbuild = .true.
                end if
                call iw_tooltip("Display the unit cell",ttshown)
             end if

             str2 = "Van der Waals Contacts [Atoms]" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                id = w%sc%get_new_representation_id()
                call w%sc%rep(id)%init(w%sc,w%view_selected,id,reptype_atoms,w%sc%style,&
                   repflavor_atoms_vdwcontacts)
                chbuild = .true.
             end if
             call iw_tooltip("Display intermolecular close contacts using van der Waals radii",ttshown)

             str2 = "Hydrogen Bonds [Atoms]" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                id = w%sc%get_new_representation_id()
                call w%sc%rep(id)%init(w%sc,w%view_selected,id,reptype_atoms,w%sc%style,&
                   repflavor_atoms_hbonds)
                chbuild = .true.
             end if
             call iw_tooltip("Display intermolecular close contacts using van der Waals radii",ttshown)

             call igEndPopup()
          end if
          call iw_tooltip("Add a new object to the view",ttshown)

          ! rest of the table
          str2 = "Objects##0,0" // c_null_char
          flags = ImGuiTableFlags_NoSavedSettings
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_NoBordersInBody)
          sz0%x = 0
          nrep = count(w%sc%rep(1:w%sc%nrep)%isinit)
          nrep = min(nrep,10)
          sz0%y = iw_calcheight(nrep,0,.true.)
          if (igBeginTable(c_loc(str2),4,flags,sz0,0._c_float)) then
             str3 = "[close button]##1closebutton" // c_null_char
             flags = ImGuiTableColumnFlags_None
             width = max(4._c_float, fontsize%y + 2._c_float)
             call igTableSetupColumn(c_loc(str3),flags,width,ic_closebutton)

             str3 = "[view button]##1viewbutton" // c_null_char
             flags = ImGuiTableColumnFlags_None
             call igTableSetupColumn(c_loc(str3),flags,0.0_c_float,ic_viewbutton)

             str3 = "[name]##1name" // c_null_char
             flags = ImGuiTableColumnFlags_WidthStretch
             call igTableSetupColumn(c_loc(str3),flags,0.0_c_float,ic_name)

             str3 = "[edit button]##1editbutton" // c_null_char
             flags = ImGuiTableColumnFlags_None
             width = iw_calcwidth(4,1)
             call igTableSetupColumn(c_loc(str3),flags,width,ic_editbutton)

             if (w%sc%representation_menu(w%id)) chbuild = .true.

             call igEndTable()
          end if

          call igEndPopup()
       end if
    end if
    call iw_tooltip("Add, remove, and modify objects",ttshown)

    ! update the draw lists and render
    if (associated(w%sc)) then
       if (chbuild .or. w%sc%timelastbuild < sysc(w%view_selected)%timelastchange) w%sc%forcebuildlists = .true.
       if (chrender .or. w%sc%forcebuildlists .or. w%sc%timelastrender < sysc(w%view_selected)%timelastchange) &
          w%forcerender = .true.

       ! continuous render if animation is active
       if (w%sc%ifreq_selected > 0.and.w%sc%iqpt_selected > 0.and.allocated(sys(w%view_selected)%c%vib).and.&
          w%sc%animation > 0) &
          w%forcerender = .true.
    end if

    ! export image
    ldum = iw_button("Tools",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft,&
       disabled=.not.associated(w%sc))
    call iw_tooltip("Show various tools operating on the view of this system",ttshown)
    if (ok) then
       str2 = "Export to Image..." // c_null_char
       enabled = associated(w%sc)
       if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,enabled)) then
          w%idexportwin = stack_create_window(wintype_exportimage,.true.,idcaller=w%id,&
             orraise=w%idexportwin)
       end if
       call iw_tooltip("Export the current view to an image file (png)",ttshown)

       str2 = "Vibrations..." // c_null_char
       if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,enabled)) then
          w%idvibrationswin = stack_create_window(wintype_vibrations,.true.,idcaller=w%id,&
             orraise=w%idvibrationswin)
       end if
       call iw_tooltip("Display an animation showing the atomic vibrations for this system",ttshown)
       call igEndPopup()
    end if


    ! camera lock
    if (w%ismain) then
       ldum = iw_button("Cam-Lock",disabled=.not.associated(w%sc),sameline=.true.,&
          popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
       if (ok) then
          str2 = "Lock All" // c_null_char
          if (igMenuItem_Bool(c_loc(str2),c_null_ptr,logical(lockbehavior==2,c_bool),.true._c_bool)) then
             lockbehavior = 2
             do k = 1, nsys
                sysc(k)%sc%lockedcam = 1
             end do
          end if
          call iw_tooltip("Lock the camera position for all loaded systems",ttshown)

          str2 = "Lock SCF Iterations Only" // c_null_char
          if (igMenuItem_Bool(c_loc(str2),c_null_ptr,logical(lockbehavior==1,c_bool),.true._c_bool)) then
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

          str2 = "Unlock All" // c_null_char
          if (igMenuItem_Bool(c_loc(str2),c_null_ptr,logical(lockbehavior==0,c_bool),.true._c_bool)) then
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
       do i = 1, win(idum)%sc%nrep
          win(idum)%sc%rep(i)%idwin = 0
       end do
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
    if (szavail%x > szavail%y) then
       ratio = szavail%x / max(szavail%y,1._c_float)
    else
       ratio = szavail%y / max(szavail%x,1._c_float)
    end if
    if (associated(w%sc)) w%sc%camratio = min(ratio,2.5_c_float)

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

    ! get hover, image rectangle coordinates, and atom idx
    idx = 0
    hover = goodsys .and. w%ilock == ilock_no
    if (hover) hover = igIsItemHovered(ImGuiHoveredFlags_None)
    call igGetItemRectMin(w%v_rmin)
    call igGetItemRectMax(w%v_rmax)
    if (hover) then
       call igGetMousePos(pos)
       call w%mousepos_to_texpos(pos)
       call w%getpixel(pos,depth,rgba)

       ! transform to atom cell ID and lattice vector
       idx(1:4) = transfer(rgba,idx(1:4))
       idx(5) = idx(1)
       if (associated(w%sc) .and. idx(1) > 0) then
          idx(1:4) = w%sc%drawlist_sph(idx(1))%idx
       else
          idx = 0
       end if
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
    if (hover .and. idx(1) > 0) then
       call igSameLine(0._c_float,-1._c_float)
       icel = idx(1)
       is = sys(w%view_selected)%c%atcel(icel)%is
       ineq = sys(w%view_selected)%c%atcel(icel)%idx
       ismol = sys(w%view_selected)%c%ismolecule

       msg = trim(sys(w%view_selected)%c%spc(is)%name)
       if (.not.ismol) then
          x0 = sys(w%view_selected)%c%atcel(icel)%x

          msg = trim(msg) // " [cellid=" // string(icel) // "+(" // string(idx(2)) // "," // string(idx(3)) //&
             "," // string(idx(4)) // "),nneqid=" // string(ineq) // ",wyckoff=" // &
             string(sys(w%view_selected)%c%at(ineq)%mult) // string(sys(w%view_selected)%c%at(ineq)%wyc)
          if (sys(w%view_selected)%c%nmol > 1) &
             msg = msg // ",molid=" // string(sys(w%view_selected)%c%idatcelmol(1,icel))
          msg = msg // "] " //&
             string(x0(1)+idx(2),'f',decimal=4) //" "// string(x0(2)+idx(3),'f',decimal=4) //" "//&
             string(x0(3)+idx(4),'f',decimal=4) // " (frac)"
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
       call w%draw_selection_tooltip(idx)

    ! Process mouse events
    call w%process_events_view(hover,idx)

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
    use gui_main, only: nsys, sysc, sys_init
    class(window), intent(inout), target :: w
    integer, intent(in) :: isys

    if (isys < 1 .or. isys > nsys) return
    if (sysc(isys)%status /= sys_init) w%forcerender = .true. ! for removing the last system in tree
    if (w%view_selected == isys) return
    w%view_selected = isys
    if (w%ismain) then
       w%sc => sysc(w%view_selected)%sc

       ! if the camera is locked, copy the camera parameters from the member
       ! of the locking group who was rendered last
       if (sysc(isys)%sc%lockedcam > 0) &
          call sysc(isys)%sc%copy_cam(idx=sysc(isys)%sc%lockedcam)
    else
       if (.not.associated(w%sc)) allocate(w%sc)
       call w%sc%end()
       call w%sc%init(w%view_selected)
    end if
    w%forcerender = .true.

  end subroutine select_view

  !> Process the mouse events in the view window
  module subroutine process_events_view(w,hover,idx)
    use interfaces_cimgui
    use scenes, only: scene, min_zoom, max_zoom
    use utils, only: translate, rotate, mult, invmult
    use tools_math, only: cross_cfloat, matinv_cfloat
    use keybindings, only: is_bind_event, is_bind_mousescroll, BIND_NAV_ROTATE,&
       BIND_NAV_ROTATE_PERP,&
       BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET, BIND_NAV_MEASURE,&
       BIND_CLOSE_FOCUSED_DIALOG
    use gui_main, only: io, nsys
    class(window), intent(inout), target :: w
    logical, intent(in) :: hover
    integer(c_int), intent(in) :: idx(5)

    type(ImVec2) :: texpos, mousepos
    real(c_float) :: ratio, pos3(3), vnew(3), vold(3), axis(3), lax
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

          call mult(xc,w%sc%world,w%sc%scenecenter)
          pos3 = w%sc%campos - xc
          pos3 = pos3 - ratio * pos3
          if (norm2(pos3) < min_zoom) &
             pos3 = pos3 / norm2(pos3) * min_zoom
          if (norm2(pos3) > max_zoom * w%sc%scenerad) &
             pos3 = pos3 / norm2(pos3) * (max_zoom * w%sc%scenerad)
          w%sc%campos = xc + pos3

          call w%sc%update_view_matrix()
          call w%sc%update_projection_matrix()
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
                w%sc%campos = xc
                call w%sc%update_view_matrix()
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
                vnew = (/texpos%x,texpos%y,w%mpos0_l(3)/)
                call w%texpos_to_view(vnew)
                pos3 = (/0._c_float,0._c_float,1._c_float/)
                axis = cross_cfloat(pos3,vnew - w%cpos0_l)
                lax = norm2(axis)
                if (lax > 1e-10_c_float) then
                   axis = axis / lax
                   call invmult(axis,w%sc%world,notrans=.true.)
                   mpos2(1) = texpos%x - w%mpos0_l(1)
                   mpos2(2) = texpos%y - w%mpos0_l(2)
                   ang = 2._c_float * norm2(mpos2) * mousesens_rot0 / w%FBOside
                   call translate(w%sc%world,w%sc%scenecenter)
                   call rotate(w%sc%world,ang,axis)
                   call translate(w%sc%world,-w%sc%scenecenter)

                   w%forcerender = .true.
                end if
                w%mpos0_l = (/texpos%x, texpos%y, 0._c_float/)
                w%cpos0_l = w%mpos0_l
                call w%texpos_to_view(w%cpos0_l)
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
                call invmult(axis,w%sc%world,notrans=.true.)
                call translate(w%sc%world,w%sc%scenecenter)
                call rotate(w%sc%world,ang,axis)
                call translate(w%sc%world,-w%sc%scenecenter)
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
    use gui_main, only: sys, fontsize
    use tools_io, only: string
    use tools_math, only: cross
    use param, only: bohrtoa, pi
    class(window), intent(inout), target :: w
    integer(c_int), intent(in) :: idx(5)

    integer :: nmsel
    integer :: msel(5,4)
    integer :: idx1(4), idx2(4), idx3(4), idx4(4)
    real*8 :: x0(3), x1(3), x2(3), d, d1, d2, ang, n0(3), n1(3)

    real(c_float), parameter :: rgbsel(4,4) = reshape((/&
       1._c_float,  0.4_c_float, 0.4_c_float, 1._c_float,&
       0.4_c_float, 1._c_float,  0.4_c_float, 1._c_float,&
       0.4_c_float, 0.4_c_float, 1._c_float, 1._c_float,&
       0.9_c_float, 0.7_c_float, 0.4_c_float, 1._c_float/),shape(rgbsel))

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
       call iw_text("1",rgba=rgbsel(:,1),sameline_nospace=.true.)
       if (nmsel > 1) then
          call iw_text("2",rgba=rgbsel(:,2),sameline_nospace=.true.)
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
          call iw_text("2",rgba=rgbsel(:,2),sameline_nospace=.true.)
          if (nmsel > 2) then
             call iw_text("3",rgba=rgbsel(:,3),sameline_nospace=.true.)
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
          call iw_text("1",rgba=rgbsel(:,1),sameline_nospace=.true.)
          call iw_text("2",rgba=rgbsel(:,2),sameline_nospace=.true.)
          if (nmsel > 2) then
             call iw_text("3",rgba=rgbsel(:,3),sameline_nospace=.true.)
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
          call iw_text("3",rgba=rgbsel(:,3),sameline_nospace=.true.)
          if (nmsel > 3) then
             call iw_text("4",rgba=rgbsel(:,4),sameline_nospace=.true.)
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
          call iw_text("2",rgba=rgbsel(:,2),sameline_nospace=.true.)
          call iw_text("3",rgba=rgbsel(:,3),sameline_nospace=.true.)
          if (nmsel > 3) then
             call iw_text("4",rgba=rgbsel(:,4),sameline_nospace=.true.)
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
       call iw_text("1",rgba=rgbsel(:,1),sameline_nospace=.true.)
       call iw_text("2",rgba=rgbsel(:,2),sameline_nospace=.true.)
       call iw_text("3",rgba=rgbsel(:,3),sameline_nospace=.true.)
       if (nmsel > 3) then
          call iw_text("4",rgba=rgbsel(:,4),sameline_nospace=.true.)
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

  !xx! edit representation

  !> Update the isys and irep in the edit represenatation window.
  module subroutine update_editrep(w)
    use windows, only: nwin, win, wintype_view
    use gui_main, only: nsys, sysc, sys_init
    class(window), intent(inout), target :: w

    integer :: isys
    logical :: doquit

    ! check the system and representation are still active
    isys = w%isys
    doquit = (isys < 1 .or. isys > nsys)
    if (.not.doquit) doquit = (sysc(isys)%status /= sys_init)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)
    if (.not.doquit) doquit = .not.(w%idparent > 0 .and. w%idparent <= nwin)
    if (.not.doquit) doquit = .not.(win(w%idparent)%isinit)
    if (.not.doquit) doquit = win(w%idparent)%type /= wintype_view

    ! if they aren't, quit the window
    if (doquit) call w%end()

  end subroutine update_editrep

  !> Draw the edit represenatation window.
  module subroutine draw_editrep(w)
    use scenes, only: representation, reptype_atoms, reptype_unitcell
    use windows, only: nwin, win, wintype_view
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG,&
       BIND_CLOSE_ALL_DIALOGS
    use gui_main, only: nsys, sysc, sys_init, g
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_calcheight, iw_checkbox
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: isys, ll, itype
    logical :: doquit, lshown, ok
    logical(c_bool) :: changed
    character(kind=c_char,len=:), allocatable, target :: str1, str2
    character(kind=c_char,len=1024), target :: txtinp
    type(ImVec2) :: szavail

    logical, save :: ttshown = .false. ! tooltip flag

    ! check the system and representation are still active
    isys = w%isys
    doquit = (isys < 1 .or. isys > nsys)
    if (.not.doquit) doquit = (sysc(isys)%status /= sys_init)
    if (.not.doquit) doquit = .not.associated(w%rep)
    if (.not.doquit) doquit = .not.w%rep%isinit
    if (.not.doquit) doquit = (w%rep%type <= 0)
    if (.not.doquit) doquit = .not.(w%idparent > 0 .and. w%idparent <= nwin)
    if (.not.doquit) doquit = .not.(win(w%idparent)%isinit)
    if (.not.doquit) doquit = .not.associated(win(w%idparent)%sc)
    if (.not.doquit) doquit = win(w%idparent)%type /= wintype_view
    if (.not.doquit) doquit = (win(w%idparent)%view_selected /= isys)

    if (.not.doquit) then
       ! whether the rep has changed
       changed = .false.

       ! system
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! name and type block
       call iw_text("Type and Name",highlight=.true.)

       ! the representation type
       itype = w%rep%type - 1
       call iw_combo_simple("##reptype","Atoms+..." // c_null_char // "Unit cell" // c_null_char,itype)
       if (w%rep%type /= itype + 1) changed = .true.
       w%rep%type = itype + 1
       call iw_tooltip("Type of object",ttshown)

       ! name text input
       str1 = "##nametextinput"
       txtinp = trim(adjustl(w%rep%name)) // c_null_char
       call igSameLine(0._c_float,-1._c_float)
       call igPushItemWidth(iw_calcwidth(30,1))
       if (igInputText(c_loc(str1),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_None,c_null_funptr,c_null_ptr)) then
          ll = index(txtinp,c_null_char)
          w%rep%name = txtinp(1:ll-1)
       end if
       call igPopItemWidth()
       call iw_tooltip("Name of this object",ttshown)

       ! shown checkbox
       changed = changed .or. iw_checkbox("Show",w%rep%shown,sameline=.true.)
       call iw_tooltip("Toggle show/hide this object",ttshown)

       ! type-dependent items
       if (w%rep%type == reptype_atoms) then
          changed = changed .or. w%draw_editrep_atoms(ttshown)
       elseif (w%rep%type == reptype_unitcell) then
          changed = changed .or. w%draw_editrep_unitcell(ttshown)
       end if

       ! rebuild draw lists if necessary
       if (changed) win(w%idparent)%sc%forcebuildlists = .true.

       ! right-align and bottom-align for the rest of the contents
       call igGetContentRegionAvail(szavail)
       call igSetCursorPosX(iw_calcwidth(7,2,from_end=.true.) - g%Style%ScrollbarSize)
       if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
          call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

       ! reset button
       if (iw_button("Reset",danger=.true.)) then
          str2 = w%rep%name
          itype = w%rep%type
          lshown = w%rep%shown
          call w%rep%init(sysc(w%isys)%sc,w%rep%id,w%rep%idrep,itype,win(w%idparent)%sc%style,&
             w%rep%flavor)
          w%rep%name = str2
          w%rep%shown = lshown
          win(w%idparent)%sc%forcebuildlists = .true.
       end if

       ! close button
       ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
       ok = ok .or. iw_button("OK",sameline=.true.)
       if (ok) doquit = .true.
    end if

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_editrep

  !> Draw the editrep (Object) window, atoms class. Returns true if
  !> the scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_atoms(w,ttshown) result(changed)
    use scenes, only: representation
    use gui_main, only: sys, g
    use tools_io, only: string
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_radiobutton, iw_calcheight, iw_clamp_color3, iw_checkbox, iw_coloredit3
    use param, only: atmcov, atmvdw, jmlcol, jmlcol2, newline
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical(c_bool) :: changed

    integer :: ispc, isys, iz, ll, ipad
    character(kind=c_char,len=1024), target :: txtinp
    character(kind=c_char,len=33), target :: txtinp2
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3, suffix
    real*8 :: x0(3)
    logical(c_bool) :: ch, ldum
    integer(c_int) :: nc(3), lst, flags, nspcpair
    real(c_float) :: sqw
    integer :: i, j, k, intable, nrow, is, ncol
    type(ImVec2) :: sz
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    integer, allocatable :: indi(:), indj(:)

    integer(c_int), parameter :: lsttrans(0:7) = (/0,1,2,2,2,3,4,5/)
    integer(c_int), parameter :: lsttransi(0:5) = (/0,1,2,5,6,7/)

    integer(c_int), parameter :: ic_sp1 = 0
    integer(c_int), parameter :: ic_sp2 = 1
    integer(c_int), parameter :: ic_shown = 2

    ! initialize
    changed = .false.
    isys = w%isys

    ! update representation to respond to changes in number of atoms and molecules
    call w%rep%update()

    ! row of display options
    changed = changed .or. iw_checkbox("Atoms##atomsglobaldisplay",w%rep%atoms_display,highlight=.true.)
    call iw_tooltip("Display atoms in the scene",ttshown)

    changed = changed .or. iw_checkbox("Bonds##bondsglobaldisplay",w%rep%bonds_display,sameline=.true.,highlight=.true.)
    call iw_tooltip("Display bonds in the scene",ttshown)

    changed = changed .or. iw_checkbox("Labels##labelsglobaldisplay",w%rep%labels_display,sameline=.true.,highlight=.true.)
    call iw_tooltip("Display atomic labels in the scene",ttshown)

    str1 = "##editrepatomstabbar" // string(w%isys) // c_null_char
    flags = ImGuiTabBarFlags_Reorderable
    flags = ior(flags,ImGuiTabBarFlags_AutoSelectNewTabs)
    if (igBeginTabBar(c_loc(str1),flags)) then

       str1 = "Atoms##editrepatoms_atomstab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
          !! atoms display !!

          ! filter
          call igAlignTextToFramePadding()
          call iw_text("Filter",highlight=.true.)
          call iw_text("(?)",sameline=.true.)
          call iw_tooltip("Show the atom if the filter expression evaluates to non-zero (true) at the atomic position. &
             &Structural variables are very useful for filters. Examples:"//newline//&
             "- '@x < 3' = all atoms with x lower than 3"//newline//&
             "- 'log($0) > 1' = log of the promolecular density higher than 1"//newline//&
             "- 'abs(@x) < 2 && abs(@y) < 2 && abs(@z) < 2' = atoms in the (-2,2) box"//newline//&
             "Click on the Help button for more info.")
          if (iw_button("Help##helpfilter",sameline=.true.)) then
             str3 = "https://aoterodelaroza.github.io/critic2/manual/arithmetics" // c_null_char
             call openLink(c_loc(str3))
          end if
          call iw_tooltip("Open the manual page about arithmetic expressions."&
             &"The 'basic usage' and 'structural variables' sections are relevant.",ttshown)

          ! filter text input
          str1 = "##filtertext" // c_null_char
          txtinp = trim(adjustl(w%rep%filter)) // c_null_char
          if (igInputText(c_loc(str1),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_EnterReturnsTrue,&
             c_null_funptr,c_null_ptr)) then
             ll = index(txtinp,c_null_char)
             w%rep%filter = txtinp(1:ll-1)

             ! test the filter
             if (sys(isys)%c%ncel > 0) then
                x0 = sys(isys)%c%atcel(1)%r
             else
                x0 = 0d0
             end if
             changed = .true.
             w%rep%errfilter = ""
          end if
          if (len_trim(w%rep%filter) == 0) w%rep%errfilter = ""
          call iw_tooltip("Apply this filter to the atoms in the system. Atoms are represented if non-zero.",&
             ttshown)
          if (iw_button("Clear",sameline=.true.)) then
             w%rep%filter = ""
             w%rep%errfilter = ""
             changed = .true.
          end if
          call iw_tooltip("Clear the filter",ttshown)
          if (len_trim(w%rep%errfilter) > 0) &
             call iw_text("Error: " // trim(w%rep%errfilter),danger=.true.)

          ! periodicity
          if (.not.sys(isys)%c%ismolecule) then
             call igAlignTextToFramePadding()
             call iw_text("Periodicity",highlight=.true.)

             ! radio buttons for the periodicity type
             changed = changed .or. iw_radiobutton("None",int=w%rep%pertype,intval=0_c_int,sameline=.true.)
             call iw_tooltip("This object is represented only in the main cell and not repeated by translation",ttshown)
             changed = changed .or. iw_radiobutton("Automatic",int=w%rep%pertype,intval=1_c_int,sameline=.true.)
             call iw_tooltip("Number of periodic cells controlled by the +/- options &
                &in the 'Scene' button of the view window",ttshown)
             changed = changed .or. iw_radiobutton("Manual",int=w%rep%pertype,intval=2_c_int,sameline=.true.)
             call iw_tooltip("Manually set the number of periodic cells",ttshown)

             ! number of periodic cells, if manual
             if (w%rep%pertype == 2_c_int) then
                ! calculate widths
                ipad = ceiling(log10(max(maxval(w%rep%ncell),1) + 0.1))
                sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
                call igPushItemWidth(sqw)

                nc = w%rep%ncell
                call igAlignTextToFramePadding()
                call iw_text("a:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##aaxis")) w%rep%ncell(1) = max(w%rep%ncell(1)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                str2 = "##aaxis" // c_null_char
                call igPushItemWidth(iw_calcwidth(ipad,1))
                ldum = igInputInt(c_loc(str2),w%rep%ncell(1),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
                call igPopItemWidth()
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##aaxis")) w%rep%ncell(1) = w%rep%ncell(1)+1

                call igSameLine(0._c_float,-1._c_float)
                call iw_text("b:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##baxis")) w%rep%ncell(2) = max(w%rep%ncell(2)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                str2 = "##baxis" // c_null_char
                call igPushItemWidth(iw_calcwidth(ipad,1))
                ldum = igInputInt(c_loc(str2),w%rep%ncell(2),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
                call igPopItemWidth()
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##baxis")) w%rep%ncell(2) = w%rep%ncell(2)+1

                call igSameLine(0._c_float,-1._c_float)
                call iw_text("c:")
                call igSameLine(0._c_float,0._c_float)
                if (iw_button("-##caxis")) w%rep%ncell(3) = max(w%rep%ncell(3)-1,1)
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                str2 = "##caxis" // c_null_char
                call igPushItemWidth(iw_calcwidth(ipad,1))
                ldum = igInputInt(c_loc(str2),w%rep%ncell(3),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
                call igPopItemWidth()
                call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
                if (iw_button("+##caxis")) w%rep%ncell(3) = w%rep%ncell(3)+1
                w%rep%ncell = max(w%rep%ncell,1)
                if (any(nc /= w%rep%ncell)) changed = .true.
                call igPopItemWidth()

                if (iw_button("Reset",sameline=.true.)) then
                   w%rep%ncell = 1
                   changed = .true.
                end if
             end if

             ! checkbox for molecular motif
             changed = changed .or. iw_checkbox("Show connected molecules",w%rep%onemotif)
             call iw_tooltip("Translate atoms to display whole molecules",ttshown)

             ! checkbox for border
             changed = changed .or. iw_checkbox("Show atoms at cell edges",w%rep%border,sameline=.true.)
             call iw_tooltip("Display atoms near the unit cell edges",ttshown)
          end if

          ! origin of the atoms
          call igPushItemWidth(iw_calcwidth(21,3))
          if (sys(isys)%c%ismolecule) then
             str1 = "Translate Origin (Å)##originatom" // c_null_char
             str2 = "%.4f" // c_null_char
             changed = changed .or. igDragFloat3(c_loc(str1),w%rep%origin,&
                0.01_c_float,-FLT_MAX,FLT_MAX,c_loc(str2),ImGuiSliderFlags_None)
             call iw_tooltip("Translation vector for the contents of the unit cell.",ttshown)
          else
             ! origin translation
             str1 = "Translate Origin (fractional)##originatom" // c_null_char
             str2 = "%.5f" // c_null_char
             changed = changed .or. igDragFloat3(c_loc(str1),w%rep%origin,&
                0.001_c_float,-FLT_MAX,FLT_MAX,c_loc(str2),ImGuiSliderFlags_None)
             call iw_tooltip("Translation vector for the contents of the unit cell.",ttshown)

             ! origin shift
             str1 = "Cell Origin Shift (fractional)##origincell" // c_null_char
             str2 = "%.5f" // c_null_char
             changed = changed .or. igDragFloat3(c_loc(str1),w%rep%tshift,&
                0.001_c_float,-FLT_MAX,FLT_MAX,c_loc(str2),ImGuiSliderFlags_None)
             call iw_tooltip("Displace the origin of the cell being represented.",ttshown)
          end if
          call igPopItemWidth()

          ! global options for atoms
          call igAlignTextToFramePadding()
          call iw_text("Global Options",highlight=.true.)
          if (iw_button("Reset##resetglobalatoms",sameline=.true.,danger=.true.)) then
             call w%rep%reset_atom_style()
             changed = .true.
          end if
          call iw_tooltip("Reset to the default settings for the atom representation")
          call iw_combo_simple("Radii ##atomradiicombo","Covalent" // c_null_char // "Van der Waals" // c_null_char,&
             w%rep%atom_radii_reset_type,changed=ch)
          call iw_tooltip("Set atomic radii to the tabulated values of this type",ttshown)

          call igSameLine(0._c_float,-1._c_float)
          call igPushItemWidth(iw_calcwidth(5,1))
          str2 = "Radii Scale##atomradiiscale" // c_null_char
          str3 = "%.3f" // c_null_char
          ch = ch .or. igDragFloat(c_loc(str2),w%rep%atom_radii_reset_scale,0.01_c_float,0._c_float,5._c_float,c_loc(str3),&
             ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Scale factor for the tabulated atomic radii",ttshown)
          call igPopItemWidth()
          if (ch) then
             do i = 1, w%rep%atom_style%ntype
                if (w%rep%atom_style%type == 0) then ! species
                   ispc = i
                elseif (w%rep%atom_style%type == 1) then ! nneq
                   ispc = sys(isys)%c%at(i)%is
                else ! ncel
                   ispc = sys(isys)%c%atcel(i)%is
                end if
                iz = sys(isys)%c%spc(ispc)%z
                if (w%rep%atom_radii_reset_type == 0) then
                   w%rep%atom_style%rad(i) = real(atmcov(iz),c_float)
                else
                   w%rep%atom_style%rad(i) = real(atmvdw(iz),c_float)
                end if
                w%rep%atom_style%rad(i) = w%rep%atom_style%rad(i) * w%rep%atom_radii_reset_scale
             end do
             changed = .true.
          end if

          ! style buttons: set color
          call iw_combo_simple("Colors ##atomcolorselect","jmol (light)" // c_null_char // "jmol2 (dark)" // c_null_char,&
             w%rep%atom_color_reset_type,changed=ch)
          call iw_tooltip("Set the color of all atoms to the tabulated values",ttshown)
          if (ch) then
             do i = 1, w%rep%atom_style%ntype
                if (w%rep%atom_style%type == 0) then ! species
                   ispc = i
                elseif (w%rep%atom_style%type == 1) then ! nneq
                   ispc = sys(isys)%c%at(i)%is
                else ! ncel
                   ispc = sys(isys)%c%atcel(i)%is
                end if
                iz = sys(isys)%c%spc(ispc)%z
                if (w%rep%atom_color_reset_type == 0) then
                   w%rep%atom_style%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
                else
                   w%rep%atom_style%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
                end if
             end do
             changed = .true.
          end if

          ! draw the atom selection widget
          changed = changed .or. atom_selection_widget(sys(isys)%c,w%rep)

          call igEndTabItem()
       end if ! begin tab item (atoms)

       str1 = "Bonds##editrepatoms_bondstab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
          !! bonds display !!
          call iw_text("(Note: the Atoms tab also controls which bonds are shown)")

          !! global options !!
          call igAlignTextToFramePadding()
          call iw_text("Global Options",highlight=.true.)
          if (iw_button("Reset##resetglobal",sameline=.true.,danger=.true.)) then
             call w%rep%reset_bond_style()
             changed = .true.
          end if
          call iw_tooltip("Reset to the covalent bonding for this system and the default settings")

          ! rest of the options (record changes)
          ch = .false.
          call igAlignTextToFramePadding()
          call iw_text("Style")
          call iw_combo_simple("##tablebondstyleglobalselect",&
             "Single color"//c_null_char//"Two colors"//c_null_char,w%rep%bond_style%style_g,sameline=.true.,changed=ch)
          call iw_tooltip("Use a single color for the bond, or two colors from the bonded atoms",ttshown)

          call iw_text(" Radius",sameline=.true.)
          str2 = "##radiusbondtableglobal" // c_null_char
          str3 = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(5,1))
          call igSameLine(0._c_float,-1._c_float)
          ch = ch .or. igDragFloat(c_loc(str2),w%rep%bond_style%rad_g,0.005_c_float,0._c_float,2._c_float,&
             c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
          call igPopItemWidth()
          call iw_tooltip("Radius of the bonds",ttshown)

          call iw_text(" Border",sameline=.true.)
          str2 = "##borderbondtableglobal" // c_null_char
          str3 = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(5,1))
          call igSameLine(0._c_float,-1._c_float)
          ch = ch .or. igDragFloat(c_loc(str2),w%rep%bond_style%border_g,0.002_c_float,0._c_float,2._c_float,&
             c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
          call igPopItemWidth()
          call iw_tooltip("Size of the border at the bond edges",ttshown)

          ! color
          call igAlignTextToFramePadding()
          call iw_text("Color")
          ch = ch .or. iw_coloredit3("##colorbondtableglobal",w%rep%bond_style%rgb_g,sameline=.true.)
          call iw_tooltip("Color of the bonds",ttshown)

          ! order
          call iw_text(" Order",sameline=.true.)
          call iw_combo_simple("##tablebondorderselectglobal",&
             "Dashed"//c_null_char//"Single"//c_null_char//"Double"//c_null_char//"Triple"//c_null_char,&
             w%rep%bond_style%order_g,sameline=.true.,changed=ldum)
          ch = ch .or. ldum
          call iw_tooltip("Bond order (dashed, single, double, etc.)",ttshown)

          ! both atoms
          call iw_text(" Both Atoms",sameline=.true.)
          ch = ch .or. iw_checkbox("##bothatomstableglobal",w%rep%bond_style%bothends_g,sameline=.true.)
          call iw_tooltip("Represent a bond if both end-atoms are in the scene (checked) or if only &
             &one end-atom is in the scene (unchecked)",ttshown)

          !! distance block !!
          call igAlignTextToFramePadding()
          call iw_text("Distances",highlight=.true.)
          call iw_text(" (",highlight=.true.,sameline_nospace=.true.)
          call iw_combo_simple("##tablebondglobaldistcombo","Factor"//c_null_char//"Range"//c_null_char,&
             w%rep%bond_style%distancetype_g,sameline_nospace=.true.)
          call iw_tooltip("Draw bonds whose lengths are a factor of the sum of atomic&
             & radii (Factor) or give bond distance range (Range)",ttshown)
          call iw_text(")",highlight=.true.,sameline_nospace=.true.)
          if (iw_button("Apply##applyglobal",sameline=.true.,danger=.true.)) then
             call w%rep%bond_style%generate_neighstars_from_globals(isys)
             w%rep%bond_style%isdef = .false.
             changed = .true.
          end if
          call iw_tooltip("Recalculate and draw bonds using the selected distance criteria",ttshown)

          if (w%rep%bond_style%distancetype_g == 0) then
             ! factor
             call igAlignTextToFramePadding()
             call iw_text("Between")
             call igSameLine(0._c_float,-1._c_float)
             call igPushItemWidth(iw_calcwidth(5,1))
             str2 = "##bondtableglobalbfmin" // c_null_char
             if (igDragFloat(c_loc(str2),w%rep%bond_style%bfmin_g,0.01_c_float,0.0_c_float,9.999_c_float,&
                c_loc(str3),ImGuiSliderFlags_AlwaysClamp)) then
             end if
             call igPopItemWidth()
             call iw_tooltip("Bonds with length below this factor times the radii are not shown",ttshown)

             call iw_text("times",sameline=.true.)
             call iw_combo_simple("##bondtableglobalradtypemin","cov."//c_null_char//"vdw"//c_null_char,&
                w%rep%bond_style%radtype_g(1),sameline=.true.)
             call iw_tooltip("Choose the atomic radii (covalent or van der Waals)",ttshown)
             call iw_text("radii",sameline=.true.)

             call igAlignTextToFramePadding()
             call iw_text("... and")
             call igSameLine(0._c_float,-1._c_float)
             call igPushItemWidth(iw_calcwidth(5,1))
             str2 = "##bondtableglobalbfmax" // c_null_char
             if (igDragFloat(c_loc(str2),w%rep%bond_style%bfmax_g,0.01_c_float,0.0_c_float,9.999_c_float,&
                c_loc(str3),ImGuiSliderFlags_AlwaysClamp)) then
             end if
             call igPopItemWidth()
             call iw_tooltip("Bonds with length above this factor times the radii are not shown",ttshown)

             call iw_text("times",sameline=.true.)
             call iw_combo_simple("##bondtableglobalradtypemax","cov."//c_null_char//"vdw"//c_null_char,&
                w%rep%bond_style%radtype_g(2),sameline=.true.)
             call iw_tooltip("Choose the atomic radii (covalent or van der Waals)",ttshown)
             call iw_text("radii",sameline=.true.)
          else
             ! range
             call igAlignTextToFramePadding()
             call iw_text("Between")
             call igSameLine(0._c_float,-1._c_float)
             call igPushItemWidth(iw_calcwidth(5,1))
             str2 = "##bondtableglobaldmin" // c_null_char
             if (igDragFloat(c_loc(str2),w%rep%bond_style%dmin_g,0.01_c_float,0.0_c_float,9.999_c_float,&
                c_loc(str3),ImGuiSliderFlags_AlwaysClamp)) then
             end if
             call igPopItemWidth()
             call iw_text("Å",sameline=.true.)
             call iw_tooltip("Bonds with length below this factor times the radii are not shown",ttshown)

             call igAlignTextToFramePadding()
             call iw_text("... and")
             call igSameLine(0._c_float,-1._c_float)
             call igPushItemWidth(iw_calcwidth(5,1))
             str2 = "##bondtableglobaldmax" // c_null_char
             if (igDragFloat(c_loc(str2),w%rep%bond_style%dmax_g,0.01_c_float,0.0_c_float,9.999_c_float,&
                c_loc(str3),ImGuiSliderFlags_AlwaysClamp)) then
             end if
             call igPopItemWidth()
             call iw_tooltip("Bonds with length above this factor times the radii are not shown",ttshown)
             call iw_text("Å",sameline=.true.)
          end if
          call igAlignTextToFramePadding()
          call iw_text("Intra/Inter-molecular")
          call iw_combo_simple("##tablebondimolselectglobal",&
             "any"//c_null_char//"intra"//c_null_char//"inter"//c_null_char,&
             w%rep%bond_style%imol_g,sameline=.true.)
          call iw_tooltip("Draw any bonds (any), only intramolecular (intra), or only intermolecular (inter)",&
             ttshown)

          !! atom selection block !!
          call iw_text("Atom Pair Selection",highlight=.true.)

          nspcpair = min(5,sys(isys)%c%nspc*(sys(isys)%c%nspc+1)/2+1)
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          str1="##tablespeciesbonding" // c_null_char
          sz%x = iw_calcwidth(17,3)
          sz%y = iw_calcheight(nspcpair,0,.false.)
          if (igBeginTable(c_loc(str1),3,flags,sz,0._c_float)) then
             ! header setup
             str2 = "Atom 1" // c_null_char
             flags = ImGuiTableColumnFlags_WidthFixed
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_sp1)

             str2 = "Atom 2" // c_null_char
             flags = ImGuiTableColumnFlags_WidthFixed
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_sp2)

             str2 = "Show" // c_null_char
             flags = ImGuiTableColumnFlags_WidthFixed
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_shown)

             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! start the clipper
             nrow = sys(isys)%c%nspc * (sys(isys)%c%nspc + 1) / 2
             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)
             allocate(indi(nrow),indj(nrow))
             k = 0
             do i = 1, sys(isys)%c%nspc
                do j = i, sys(isys)%c%nspc
                   k = k + 1
                   indi(k) = i
                   indj(k) = j
                end do
             end do

             ! draw the rows
             do while(ImGuiListClipper_Step(clipper))
                call c_f_pointer(clipper,clipper_f)
                do k = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                   i = indi(k)
                   j = indj(k)

                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   suffix = "_" // string(i) // "_" // string(j)

                   ! species
                   if (igTableSetColumnIndex(ic_sp1)) then
                      call igAlignTextToFramePadding()
                      call iw_text(trim(sys(isys)%c%spc(i)%name))
                   end if
                   if (igTableSetColumnIndex(ic_sp2)) &
                      call iw_text(trim(sys(isys)%c%spc(j)%name))

                   ! shown
                   if (igTableSetColumnIndex(ic_shown)) then
                      if (iw_checkbox("##bondtableshown" // suffix,w%rep%bond_style%shown_g(i,j))) then
                         ch = .true.
                         w%rep%bond_style%shown_g(j,i) = w%rep%bond_style%shown_g(i,j)
                      end if
                      call iw_tooltip("Toggle display of bonds connecting these atom types",ttshown)
                   end if
                end do ! clipper range
             end do ! clipper step

             ! end the clipper and the table
             deallocate(indi,indj)
             call ImGuiListClipper_End(clipper)
             call igEndTable()
          end if ! begintable

          ! immediately update if non-distances have changed
          if (ch) changed = .true.

          call igEndTabItem()
       end if ! begin tab item (bonds)

       str1 = "Labels##editrepatoms_labelstab" // c_null_char
       flags = ImGuiTabItemFlags_None
       if (igBeginTabItem(c_loc(str1),c_null_ptr,flags)) then
          call iw_text("(Note: the Atoms tab also controls which labels are shown)")
          !! labels display !!

          ! label styles
          call iw_text("Global Options",highlight=.true.)
          if (sys(isys)%c%ismolecule) then
             lst = lsttrans(w%rep%label_style%style)
             call iw_combo_simple("Text##labelcontentselect","Atomic symbol"//c_null_char//&
                "Atom name"// c_null_char//"Atom ID"// c_null_char//&
                "Species ID"// c_null_char// "Atomic number"// c_null_char// "Molecule ID"// c_null_char,&
                lst,changed=ch)
             w%rep%label_style%style = lsttransi(lst)
          else
             call iw_combo_simple("Text##labelcontentselect","Atomic symbol"//c_null_char//&
                "Atom name"//c_null_char//"Cell atom ID"//c_null_char//&
                "Cell atom ID + lattice vector"//c_null_char//"Symmetry-unique atom ID"//c_null_char//&
                "Species ID"//c_null_char//"Atomic number"//c_null_char//"Molecule ID"//c_null_char//&
                "Wyckoff position"//c_null_char,&
                w%rep%label_style%style,changed=ch)
          end if
          if (ch) call w%rep%reset_label_style()
          call iw_tooltip("Text to display in the atom labels",ttshown)
          changed = changed .or. ch

          ! scale, constant size, color
          str2 = "Scale##labelscale" // c_null_char
          str3 = "%.2f" // c_null_char
          call igPushItemWidth(iw_calcwidth(4,1))
          changed = changed .or. igDragFloat(c_loc(str2),w%rep%label_style%scale,0.01_c_float,&
             0._c_float,10._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
          call igPopItemWidth()
          call iw_tooltip("Scale factor for the atom labels",ttshown)

          changed = changed .or. iw_checkbox("Constant size##labelconstsize",&
             w%rep%label_style%const_size,sameline=.true.)
          call iw_tooltip("Labels have constant size (on) or labels scale with the&
             & size of the associated atom (off)",ttshown)

          changed = changed .or. iw_coloredit3("Color##labelcolor",w%rep%label_style%rgb,sameline=.true.)
          call iw_tooltip("Color of the atom labels",ttshown)

          ! offset
          call igPushItemWidth(iw_calcwidth(21,3))
          str2 = "Offset (Å)" // c_null_char
          str3 = "%.3f" // c_null_char
          changed = changed .or. igDragFloat3(c_loc(str2),w%rep%label_style%offset,&
             0.001_c_float,99.999_c_float,99.999_c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
          call iw_tooltip("Offset the position of the labels relative to the atom center",ttshown)
          call igPopItemWidth()

          ! table for label selection
          call iw_text("Label Selection",highlight=.true.)

          ! number of entries in the table
          select case(w%rep%label_style%style)
          case (0,1,5,6)
             intable = 0 ! species
             nrow = sys(isys)%c%nspc
             ncol = 5
             call iw_text("(per species)",sameline=.true.)
          case (2,3)
             intable = 1 ! complete cell list
             nrow = sys(isys)%c%ncel
             ncol = 5
             call iw_text("(per atom)",sameline=.true.)
          case (4,8)
             intable = 2 ! non-equivalent list
             nrow = sys(isys)%c%nneq
             ncol = 5
             call iw_text("(per symmetry-unique atom)",sameline=.true.)
          case (7)
             intable = 3 ! molecules
             nrow = sys(isys)%c%nmol
             ncol = 3
             call iw_text("(per molecule)",sameline=.true.)
          end select

          ! the table itself
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          str1="##tablespecieslabels" // c_null_char
          sz%x = iw_calcwidth(30,ncol)
          sz%y = iw_calcheight(min(8,nrow+1),0,.false.)
          if (igBeginTable(c_loc(str1),ncol,flags,sz,0._c_float)) then
             ncol = -1

             ! header setup
             str2 = "Id" // c_null_char
             ncol = ncol + 1
             flags = ImGuiTableColumnFlags_WidthFixed
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

             if (intable < 3) then
                str2 = "Atom" // c_null_char
                ncol = ncol + 1
                flags = ImGuiTableColumnFlags_WidthFixed
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

                str2 = "Z " // c_null_char
                ncol = ncol + 1
                flags = ImGuiTableColumnFlags_WidthFixed
                call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)
             end if

             str2 = "Show" // c_null_char
             ncol = ncol + 1
             flags = ImGuiTableColumnFlags_WidthFixed
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

             str2 = "Text" // c_null_char
             ncol = ncol + 1
             flags = ImGuiTableColumnFlags_WidthStretch
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ncol)

             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! start the clipper
             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)

             ! draw the rows
             do while(ImGuiListClipper_Step(clipper))
                call c_f_pointer(clipper,clipper_f)
                do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd

                   ! set up the next row
                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   suffix = "_" // string(i)
                   ncol = -1

                   ! id
                   ncol = ncol + 1
                   if (igTableSetColumnIndex(ncol)) then
                      call igAlignTextToFramePadding()
                      call iw_text(string(i))
                   end if

                   if (intable == 0) then ! species
                      is = i
                   elseif (intable == 1) then ! complete cell list
                      is = sys(isys)%c%atcel(i)%is
                   elseif (intable == 2) then ! non-equivalent list
                      is = sys(isys)%c%at(i)%is
                   end if

                   ! atom
                   if (intable < 3) then
                      ncol = ncol + 1
                      if (igTableSetColumnIndex(ncol)) &
                         call iw_text(trim(sys(isys)%c%spc(is)%name))

                      ! Z
                      ncol = ncol + 1
                      if (igTableSetColumnIndex(ncol)) &
                         call iw_text(string(sys(isys)%c%spc(is)%z))
                   end if

                   ! shown
                   ncol = ncol + 1
                   if (igTableSetColumnIndex(ncol)) then
                      changed = changed .or. iw_checkbox("##labeltableshown" // suffix,w%rep%label_style%shown(i))
                      call iw_tooltip("Toggle display of labels for these atoms/molecules",ttshown)
                   end if

                   ! text
                   ncol = ncol + 1
                   if (igTableSetColumnIndex(ncol)) then
                      str1 = "##labeltabletext" // suffix // c_null_char
                      txtinp2 = trim(w%rep%label_style%str(i)) // c_null_char
                      call igPushItemWidth(iw_calcwidth(15,1))
                      if (igInputText(c_loc(str1),c_loc(txtinp2),32_c_size_t,ImGuiInputTextFlags_None,&
                         c_null_funptr,c_null_ptr)) then
                         ll = index(txtinp2,c_null_char)
                         w%rep%label_style%str(i) = txtinp2(1:ll-1)
                         changed = .true.
                      end if
                      call igPopItemWidth()
                      call iw_tooltip("Text for the atomic labels",ttshown)
                   end if
                end do ! table rows: clipper range
             end do ! table rows: clipper step

             ! end the clipper and the table
             call ImGuiListClipper_End(clipper)
             call igEndTable()
          end if ! begintable

          call igEndTabItem()
       end if ! begin tab item (labels)
       call igEndTabBar()
    end if ! begin tab bar

  end function draw_editrep_atoms

  !> Draw the editrep window, unit cell class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_unitcell(w,ttshown) result(changed)
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_radiobutton, iw_button,&
       iw_clamp_color3, iw_checkbox, iw_coloredit3
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical(c_bool) :: changed, ldum
    integer :: ipad
    real(c_float) :: sqw
    integer(c_int) :: nc(3)

    character(kind=c_char,len=:), allocatable, target :: str1, str2
    logical :: ch

    ! initialize
    changed = .false.

    ! periodicity
    call igAlignTextToFramePadding()
    call iw_text("Periodicity",highlight=.true.)

    ! radio buttons for the periodicity type
    changed = changed .or. iw_radiobutton("None",int=w%rep%pertype,intval=0_c_int,sameline=.true.)
    call iw_tooltip("This object is represented only in the main cell and not repeated by translation",ttshown)
    changed = changed .or. iw_radiobutton("Automatic",int=w%rep%pertype,intval=1_c_int,sameline=.true.)
    call iw_tooltip("Number of periodic cells controlled by the +/- options in the view menu",ttshown)
    changed = changed .or. iw_radiobutton("Manual",int=w%rep%pertype,intval=2_c_int,sameline=.true.)
    call iw_tooltip("Manually set the number of periodic cells",ttshown)

    ! number of periodic cells, if manual
    if (w%rep%pertype == 2_c_int) then
       ! calculate widths
       ipad = ceiling(log10(max(maxval(w%rep%ncell),1) + 0.1))
       sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
       call igPushItemWidth(sqw)

       nc = w%rep%ncell
       call igAlignTextToFramePadding()
       call iw_text("a:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##aaxis")) w%rep%ncell(1) = max(w%rep%ncell(1)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       str2 = "##aaxis" // c_null_char
       call igPushItemWidth(iw_calcwidth(ipad,1))
       ldum = igInputInt(c_loc(str2),w%rep%ncell(1),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
       call igPopItemWidth()
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##aaxis")) w%rep%ncell(1) = w%rep%ncell(1)+1

       call igSameLine(0._c_float,-1._c_float)
       call iw_text("b:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##baxis")) w%rep%ncell(2) = max(w%rep%ncell(2)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       str2 = "##baxis" // c_null_char
       call igPushItemWidth(iw_calcwidth(ipad,1))
       ldum = igInputInt(c_loc(str2),w%rep%ncell(2),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
       call igPopItemWidth()
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##baxis")) w%rep%ncell(2) = w%rep%ncell(2)+1

       call igSameLine(0._c_float,-1._c_float)
       call iw_text("c:")
       call igSameLine(0._c_float,0._c_float)
       if (iw_button("-##caxis")) w%rep%ncell(3) = max(w%rep%ncell(3)-1,1)
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       str2 = "##caxis" // c_null_char
       call igPushItemWidth(iw_calcwidth(ipad,1))
       ldum = igInputInt(c_loc(str2),w%rep%ncell(3),-1_c_int,-100_c_int,ImGuiInputTextFlags_EnterReturnsTrue)
       call igPopItemWidth()
       call igSameLine(0._c_float,0.5_c_float*g%Style%FramePadding%x)
       if (iw_button("+##caxis")) w%rep%ncell(3) = w%rep%ncell(3)+1
       w%rep%ncell = max(w%rep%ncell,1)
       if (any(nc /= w%rep%ncell)) changed = .true.
       call igPopItemWidth()

       if (iw_button("Reset",sameline=.true.)) then
          w%rep%ncell = 1
          changed = .true.
       end if
    end if

    !! styles
    call iw_text("Style",highlight=.true.)
    changed = changed .or. iw_checkbox("Color crystallographic axes",w%rep%uc_coloraxes)
    call iw_tooltip("Represent crystallographic axes with colors (a=red,b=green,c=blue)",ttshown)

    str1 = "Radius##outer" // c_null_char
    str2 = "%.3f" // c_null_char
    call igPushItemWidth(iw_calcwidth(5,1))
    changed = changed .or. igDragFloat(c_loc(str1),w%rep%uc_radius,0.005_c_float,0._c_float,&
       5._c_float,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
    w%rep%uc_radius = max(w%rep%uc_radius,0._c_float)
    call igPopItemWidth()
    call iw_tooltip("Radii of the unit cell edges",ttshown)

    ch = iw_coloredit3("Color",w%rep%uc_rgb,sameline=.true.)
    call iw_tooltip("Color of the unit cell edges",ttshown)
    if (ch) then
       w%rep%uc_rgb = min(w%rep%uc_rgb,1._c_float)
       w%rep%uc_rgb = max(w%rep%uc_rgb,0._c_float)
       changed = .true.
    end if

    !! inner divisions
    call iw_text("Inner Divisions",highlight=.true.)
    changed = changed .or. iw_checkbox("Display inner divisions",w%rep%uc_inner)
    call iw_tooltip("Represent the inner divisions inside a supercell",ttshown)
    if (w%rep%uc_inner) then
       str1 = "Radius##inner" // c_null_char
       str2 = "%.3f" // c_null_char
       call igPushItemWidth(iw_calcwidth(5,1))
       changed = changed .or. igDragFloat(c_loc(str1),w%rep%uc_radiusinner,0.005_c_float,0._c_float,&
          5._c_float,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
       w%rep%uc_radiusinner = max(w%rep%uc_radiusinner,0._c_float)
       call igPopItemWidth()
       call iw_tooltip("Radii of the inner unit cell edges",ttshown)

       changed = changed .or. iw_checkbox("Use dashed lines",w%rep%uc_innerstipple)
       call iw_tooltip("Use dashed lines for the inner cell divisions",ttshown)

       if (w%rep%uc_innerstipple) then
          str1 = "Dash length (Å)" // c_null_char
          str2 = "%.1f" // c_null_char
          call igPushItemWidth(iw_calcwidth(5,1))
          changed = changed .or. igDragFloat(c_loc(str1),w%rep%uc_innersteplen,0.1_c_float,0._c_float,&
             100._c_float,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
          w%rep%uc_innersteplen = max(w%rep%uc_innersteplen,0._c_float)
          call igPopItemWidth()
          call iw_tooltip("Length of the dashed lines for the inner cell divisions (in Å)",ttshown)
       end if
    end if

    ! origin of the unit cell
    call iw_text("Origin Shift",highlight=.true.)
    str1 = "##originucx" // c_null_char
    str2 = "%.5f" // c_null_char
    call igPushItemWidth(iw_calcwidth(21,3))
    changed = changed .or. igDragFloat3(c_loc(str1),w%rep%origin,&
       0.001_c_float,-FLT_MAX,FLT_MAX,c_loc(str2),ImGuiSliderFlags_None)
    call iw_tooltip("Coordinates for the origin shift of the unit cell",ttshown)
    call igPopItemWidth()

  end function draw_editrep_unitcell

  !> Draw the export image window
  module subroutine draw_exportimage(w)
    use interfaces_opengl3
    use interfaces_stb
    use gui_main, only: sysc, sys_init, nsys, g
    use windows, only: wintype_dialog, wpurp_dialog_saveimagefile
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip, get_current_working_dir,&
       iw_checkbox
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG,&
       BIND_CLOSE_ALL_DIALOGS
    use tools_io, only: ferror, string
    use param, only: dirsep
    class(window), intent(inout), target :: w

    logical :: doquit, ok, goodsys, okvalid
    integer :: oid, isys, width, height
    integer(c_int), target :: msFBO, endFBO ! framebuffer
    integer(c_int), target :: msFBOdepth, endFBOdepth ! framebuffer, depth buffer
    integer(c_int), target :: msFBOtex, endFBOtex ! framebuffer, texture
    integer(c_signed_char), allocatable, target :: data(:)
    integer(c_int) :: idum, origin(2)
    character(kind=c_char,len=:), allocatable, target :: str, str1, str2
    type(ImVec2) :: x0, x1, szavail
    logical(c_bool) :: ldum

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize state
    if (w%firstpass) then
       w%okfile = get_current_working_dir() // dirsep // "image.png"
       w%okfilter = "PNG"
       w%nsample = 16
       w%jpgquality = 90
       w%exportview = .true.
       w%npixel = win(w%idparent)%FBOside
       w%transparentbg = .true.
       w%errmsg = ""
    end if

    ! initialize
    doquit = .false.
    if (associated(win(w%idparent)%sc)) then
       isys = win(w%idparent)%sc%id
    else
       isys = win(w%idparent)%view_selected
       doquit = .true.
    end if

    ! check if we have info from the export image window when it
    ! closes and recover it
    call update_window_id(w%idsave,oid)
    if (oid /= 0) then
       if (win(oid)%okfile_set) then
          w%okfile = win(oid)%okfile
          w%okfilter = win(oid)%okfilter
       end if
    end if

    ! Image file button
    call iw_text("Image File",highlight=.true.)
    if (iw_button("File",disabled=(w%idsave > 0),danger=.true.)) &
       w%idsave = stack_create_window(wintype_dialog,.true.,wpurp_dialog_saveimagefile)
    call iw_tooltip("Choose the file to save the image to",ttshown)
    call iw_text(w%okfile,sameline=.true.)

    ! render settings
    call iw_text("Render Settings",highlight=.true.)
    str1 = "Render buffer size (pixels)" // c_null_char
    str2 = "%d" // c_null_char
    call igPushItemWidth(iw_calcwidth(5,1))
    ldum = igDragInt(c_loc(str1),w%npixel,10._c_float,1_c_int,81920_c_int,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
    call igPopItemWidth()
    call iw_tooltip("Size of the render buffer in pixels",ttshown)

    str1 = "Number of samples for anti-aliasing" // c_null_char
    str2 = "%d" // c_null_char
    call igPushItemWidth(iw_calcwidth(2,1))
    ldum = igDragInt(c_loc(str1),w%nsample,1._c_float,1_c_int,32_c_int,c_loc(str2),&
       ImGuiSliderFlags_AlwaysClamp)
    call igPopItemWidth()
    call iw_tooltip("Number of samples for the anti-aliasing render",ttshown)

    ldum = iw_checkbox("Export the viewport only",w%exportview)
    call iw_tooltip("Export the viewport only or the whole render buffer",ttshown)

    ldum = iw_checkbox("Transparent background",w%transparentbg)
    call iw_tooltip("Make the background transparent in the exported image",ttshown)

    ! image settings
    if (w%okfilter(1:3) == "JPE") then
       call iw_text("Image Settings",highlight=.true.)
       str1 = "JPEG Quality" // c_null_char
       str2 = "%d" // c_null_char
       call igPushItemWidth(iw_calcwidth(3,1))
       ldum = igDragInt(c_loc(str1),w%jpgquality,1._c_float,1_c_int,100_c_int,c_loc(str2),&
          ImGuiSliderFlags_AlwaysClamp)
       call igPopItemWidth()
       call iw_tooltip("Quality and weight of the JPEG file",ttshown)
    end if

    ! maybe the error message
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! final buttons: OK
    okvalid = (len_trim(w%okfile) > 0)
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) .and. okvalid
    ok = ok .or. iw_button("OK",disabled=.not.okvalid)
    if (ok) then
       ! reset the error message
       w%errmsg = ""

       ! generate textures and buffers
       call glGenTextures(1, c_loc(msFBOtex))
       call glGenRenderbuffers(1, c_loc(msFBOdepth))
       call glGenFramebuffers(1, c_loc(msFBO))
       call glGenTextures(1, c_loc(endFBOtex))
       call glGenRenderbuffers(1, c_loc(endFBOdepth))
       call glGenFramebuffers(1, c_loc(endFBO))

       ! textures
       call glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, msFBOtex)
       call glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, w%nsample, GL_RGBA, w%npixel, w%npixel,&
          int(GL_TRUE,c_signed_char))
       call glBindTexture(GL_TEXTURE_2D, 0)
       call glBindTexture(GL_TEXTURE_2D, endFBOtex)
       call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w%npixel, w%npixel, 0, GL_RGBA, GL_UNSIGNED_BYTE, c_null_ptr)
       call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
       call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
       call glBindTexture(GL_TEXTURE_2D, 0)

       ! render buffer
       call glBindRenderbuffer(GL_RENDERBUFFER, msFBOdepth)
       call glRenderbufferStorageMultisample(GL_RENDERBUFFER, w%nsample, GL_DEPTH_COMPONENT, w%npixel, w%npixel)
       call glBindRenderbuffer(GL_RENDERBUFFER, 0)
       call glBindRenderbuffer(GL_RENDERBUFFER, endFBOdepth)
       call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w%npixel, w%npixel)
       call glBindRenderbuffer(GL_RENDERBUFFER, 0)

       ! frame buffer
       call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
       call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, msFBOtex, 0)
       call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, msFBOdepth)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
          w%errmsg = "Failed creating the multi-sampled framebuffer (too large?)"
          goto 999
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, endFBO)
       call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, endFBOtex, 0)
       call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, endFBOdepth)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
          w%errmsg = "Failed creating the render framebuffer (too large?)"
          goto 999
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! render the scene to the multisampled framebuffer
       call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
       call glViewport(0_c_int,0_c_int,w%npixel,w%npixel)
       if (w%transparentbg) then
          call glClearColor(win(w%idparent)%sc%bgcolor(1),win(w%idparent)%sc%bgcolor(2),&
             win(w%idparent)%sc%bgcolor(3),0._c_float)
       else
          call glClearColor(win(w%idparent)%sc%bgcolor(1),win(w%idparent)%sc%bgcolor(2),&
             win(w%idparent)%sc%bgcolor(3),1._c_float)
       end if
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
       goodsys = (isys >= 1 .and. isys <= nsys)
       if (goodsys) goodsys = (sysc(isys)%status == sys_init)
       if (goodsys) call win(w%idparent)%sc%render()
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! blit the multisampled buffer to the normal colorbuffer
       call glBindFramebuffer(GL_READ_FRAMEBUFFER, msFBO)
       call glBindFramebuffer(GL_DRAW_FRAMEBUFFER, endFBO)
       call glBlitFramebuffer(0, 0, w%npixel, w%npixel, 0, 0, w%npixel, w%npixel, GL_COLOR_BUFFER_BIT, GL_LINEAR)
       call glBindFramebuffer(GL_READ_FRAMEBUFFER, 0)
       call glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0)

       ! Read from the regular framebuffer into the data array
       call glBindFramebuffer(GL_FRAMEBUFFER, endFBO)
       if (.not.w%exportview) then
          ! whole texture
          width = w%npixel
          height = w%npixel
          allocate(data(4 * width * height))
          data = 0_c_signed_char
          call glReadPixels(0, 0, w%npixel, w%npixel, GL_RGBA, GL_UNSIGNED_BYTE, c_loc(data))
       else
          ! viewport only
          x0%x = win(w%idparent)%v_rmin%x
          x0%y = win(w%idparent)%v_rmin%y
          x1%x = win(w%idparent)%v_rmax%x
          x1%y = win(w%idparent)%v_rmax%y
          call win(w%idparent)%mousepos_to_texpos(x0)
          call win(w%idparent)%mousepos_to_texpos(x1)
          width = min(nint((x1%x - x0%x) / real(win(w%idparent)%FBOside,8) * w%npixel),w%npixel)
          height = min(nint((x1%y - x0%y) / real(win(w%idparent)%FBOside,8) * w%npixel),w%npixel)
          origin(1) = max(nint(x0%x / real(win(w%idparent)%FBOside,8) * w%npixel),0)
          origin(2) = max(nint(x0%y / real(win(w%idparent)%FBOside,8) * w%npixel),0)
          allocate(data(4 * width * height))
          data = 0_c_signed_char
          call glReadPixels(origin(1), origin(2), width, height, GL_RGBA, GL_UNSIGNED_BYTE, c_loc(data))
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
          w%errmsg = "Error rendering export image"

       ! write the file
       str = trim(w%okfile) // c_null_char
       if (w%okfilter(1:3) == "PNG") then
          idum = stbi_write_png(c_loc(str), width, height, 4, c_loc(data), 4*width)
       elseif (w%okfilter(1:3) == "BMP") then
          idum = stbi_write_bmp(c_loc(str), width, height, 4, c_loc(data))
       elseif (w%okfilter(1:3) == "TGA") then
          idum = stbi_write_tga(c_loc(str), width, height, 4, c_loc(data))
       elseif (w%okfilter(1:3) == "JPE") then
          idum = stbi_write_jpg(c_loc(str), width, height, 4, c_loc(data), w%jpgquality)
       end if
       if (idum == 0) &
          w%errmsg = "Error exporting image to file: "  // string(w%okfile)

999    continue

       ! delete the buffers
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)
       call glDeleteTextures(1, c_loc(msFBOtex))
       call glDeleteRenderbuffers(1, c_loc(msFBOdepth))
       call glDeleteFramebuffers(1, c_loc(msFBO))
       call glDeleteTextures(1, c_loc(endFBOtex))
       call glDeleteRenderbuffers(1, c_loc(endFBOdepth))
       call glDeleteFramebuffers(1, c_loc(endFBO))

       ! quit if no error message
       if (len_trim(w%errmsg) == 0) doquit = .true.
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_exportimage

  !> Draw the vibrations window
  module subroutine draw_vibrations(w)
    use crystalseedmod, only: crystalseed
    use scenes, only: anim_speed_default, anim_amplitude_default, anim_amplitude_max,&
       anim_speed_max
    use gui_main, only: sysc, sys, nsys, sys_init, g, add_systems_from_seeds,&
       launch_initialization_thread, time
    use utils, only: iw_text, iw_button, iw_tooltip, iw_calcheight, iw_calcwidth,&
       iw_combo_simple, iw_radiobutton
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use tools_math, only: rational_approx
    use tools_io, only: string, ioj_right, uout
    use param, only: cm1tothz, bohrtoa
    class(window), intent(inout), target :: w

    logical(c_bool) :: selected
    logical :: doquit, system_ok, vib_ok, goodparent, ldum, fset
    integer :: isys, oid, i, digits, idx
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, strl
    type(ImVec2) :: sz0, szero, szavail
    real*8 :: unitfactor, xx(3)
    integer*8 :: q, r(3)
    type(crystalseed), allocatable :: seed(:)

    integer, parameter :: ic_q_id = 0
    integer, parameter :: ic_q_qpt = 1
    real*8, parameter :: rational_approx_eps = 1d-3

    logical, save :: ttshown = .false. ! tooltip flag

    ! do we have a good parent window?
    goodparent = w%idparent > 0 .and. w%idparent <= nwin
    if (goodparent) goodparent = win(w%idparent)%isinit
    if (goodparent) goodparent = (win(w%idparent)%type == wintype_view)

    ! initialize state
    if (w%firstpass) then
       w%errmsg = ""
       if (goodparent) then
          if (associated(win(w%idparent)%sc)) then
             win(w%idparent)%sc%iqpt_selected = 0
             win(w%idparent)%sc%ifreq_selected = 0
             win(w%idparent)%sc%animation = 0
          end if
       end if
       w%ifrequnit = 0
       w%iqptunit = 0
    end if

    ! initialize
    isys = 0
    szero%x = 0
    szero%y = 0
    doquit = .false.
    if (.not.doquit) doquit = .not.goodparent
    if (.not.doquit) then
       if (associated(win(w%idparent)%sc)) then
          isys = win(w%idparent)%sc%id
       else
          isys = win(w%idparent)%view_selected
          doquit = .true.
       end if
    end if
    system_ok = (isys > 0 .and. isys <= nsys)
    if (system_ok) system_ok = (sysc(isys)%status == sys_init)

    ! check if we have info from the export image window when it
    ! closes and recover it
    call update_window_id(w%idsave,oid)
    if (oid /= 0) then
       if (system_ok .and. win(oid)%okfile_set) then
          str1 = win(oid)%okfile
          do while (.true.)
             idx = index(str1,c_null_char)
             if (idx == 0) exit
             call sys(isys)%c%read_vibrations_file(str1(1:idx-1),win(oid)%dialog_data%isformat,w%errmsg)
             if (len_trim(w%errmsg) > 0) then
                write (uout,'(A)') w%errmsg
                w%errmsg = ""
             end if
             str1 = str1(idx+1:)
          end do
       end if
    end if

    ! vibrations ok?
    vib_ok = system_ok
    if (vib_ok) vib_ok = allocated(sys(isys)%c%vib)
    if (vib_ok) vib_ok = (sys(isys)%c%vib%nqpt > 0) .and. (sys(isys)%c%vib%nfreq > 0)
    if (vib_ok) vib_ok = associated(win(w%idparent)%sc)

    ! header
    if (system_ok) then
       ! system name
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! source of vibration data
       call igAlignTextToFramePadding()
       call iw_text("Vibration data",highlight=.true.)

       if (iw_button("Clear",sameline=.true.,danger=.true.)) then
          call sys(isys)%c%clear_vibrations()
          win(w%idparent)%sc%iqpt_selected = 0
          win(w%idparent)%sc%ifreq_selected = 0
          vib_ok = .false.
       end if
       call iw_tooltip("Clear vibration data for this system",ttshown)

       if (.not.vib_ok) then
          if (iw_button("Load",danger=.true.,sameline=.true.)) &
             w%idsave = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openvibfile)
          call iw_tooltip("Load vibration data from a file for this system",ttshown)
          call iw_text("<none>",sameline=.true.)
       else
          call iw_text(sys(isys)%c%vib%file,sameline=.true.)
       end if
    end if

    ! maybe the error message
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! rest of animation stuff
    if (vib_ok) then

       if (.not.sys(isys)%c%ismolecule) then
          ! q-points table
          call igBeginGroup()
          call igAlignTextToFramePadding()
          call iw_text("Q-points",highlight=.true.)
          call iw_combo_simple("##qptunit","fractional" // c_null_char // "1/bohr" // c_null_char //&
             "1/Å" // c_null_char,w%iqptunit,sameline=.true.)
          call iw_tooltip("Units for the q-points",ttshown)

          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          str1="##tablevibrationqpoints" // c_null_char
          sz0%x = iw_calcwidth(37,0)
          sz0%y = iw_calcheight(5,0,.false.)
          if (igBeginTable(c_loc(str1),2,flags,sz0,0._c_float)) then
             ! header setup
             str2 = "Id" // c_null_char
             flags = ImGuiTableColumnFlags_WidthFixed
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_q_id)

             str2 = "Coordinates" // c_null_char
             flags = ImGuiTableColumnFlags_WidthFixed
             call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_q_qpt)
             call igTableSetupScrollFreeze(0, 1) ! top row always visible

             ! draw the header
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             ! draw the rows
             do i = 1, sys(isys)%c%vib%nqpt
                call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)

                xx = sys(isys)%c%vib%qpt(:,i)
                if (w%iqptunit == 0) then ! fractional
                   digits = 5
                elseif (w%iqptunit == 1) then ! 1/bohr
                   xx = sys(isys)%c%rx2rc(xx)
                   digits = 6
                else ! 1/ang
                   xx = sys(isys)%c%rx2rc(xx) / bohrtoa
                   digits = 6
                end if

                ! id
                if (igTableSetColumnIndex(ic_q_id)) then
                   ! selectable
                   call igAlignTextToFramePadding()
                   strl = "##selectq" // string(i) // c_null_char
                   flags = ImGuiSelectableFlags_SpanAllColumns
                   flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
                   selected = (win(w%idparent)%sc%iqpt_selected == i)
                   if (igSelectable_Bool(c_loc(strl),selected,flags,szero)) then
                      win(w%idparent)%sc%iqpt_selected = i
                      win(w%idparent)%sc%forcebuildlists = .true.
                   end if

                   ! text
                   call iw_text(string(i),sameline=.true.)
                end if

                ! coordinates
                if (igTableSetColumnIndex(ic_q_qpt)) then
                   s = string(xx(1),'f',length=10,decimal=digits,justify=ioj_right)//&
                      string(xx(2),'f',length=10,decimal=digits,justify=ioj_right)//&
                      string(xx(3),'f',length=10,decimal=digits,justify=ioj_right)
                   call iw_text(s)
                end if
             end do ! i = 1, sys(isys)%c%vib%nqpt
             call igEndTable()
          end if ! begintable
          call igEndGroup()
          call igSameLine(0._c_float,-1._c_float)
       else
          win(w%idparent)%sc%iqpt_selected = 1
       end if

       ! frequency table
       call igBeginGroup()
       ! header
       call igAlignTextToFramePadding()
       call iw_text("Frequencies",highlight=.true.)
       call iw_combo_simple("##frequnit","1/cm" // c_null_char // "THz" // c_null_char,&
          w%ifrequnit,sameline=.true.)
       call iw_tooltip("Units for the frequencies",ttshown)
       if (w%ifrequnit == 0) then
          unitfactor = 1d0
          digits = 2
       else
          unitfactor = cm1tothz
          digits = 4
       end if

       ! frequencies table
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       str1="##tablevibrationfreqs" // c_null_char
       sz0%x = iw_calcwidth(20,0)
       sz0%y = iw_calcheight(5,0,.false.)
       if (igBeginTable(c_loc(str1),2,flags,sz0,0._c_float)) then
          ! header setup
          str2 = "Id" // c_null_char
          flags = ImGuiTableColumnFlags_WidthFixed
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_q_id)

          str2 = "Frequency" // c_null_char
          flags = ImGuiTableColumnFlags_WidthFixed
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_q_qpt)
          call igTableSetupScrollFreeze(0, 1) ! top row always visible

          ! draw the header
          call igTableHeadersRow()
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())

          ! check if the qpt/frequency has been set
          fset = (win(w%idparent)%sc%iqpt_selected > 0 .and. win(w%idparent)%sc%ifreq_selected > 0)

          if (win(w%idparent)%sc%iqpt_selected > 0) then
             ! draw the rows
             do i = 1, sys(isys)%c%vib%nfreq
                call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)

                ! id
                if (igTableSetColumnIndex(ic_q_id)) then
                   ! selectable
                   call igAlignTextToFramePadding()
                   strl = "##selectf" // string(i) // c_null_char
                   flags = ImGuiSelectableFlags_SpanAllColumns
                   flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
                   selected = (win(w%idparent)%sc%ifreq_selected == i)
                   if (igSelectable_Bool(c_loc(strl),selected,flags,szero)) then
                      win(w%idparent)%sc%ifreq_selected = i
                      win(w%idparent)%sc%forcebuildlists = .true.
                   end if

                   ! text
                   call iw_text(string(i),sameline=.true.)
                end if

                ! frequency
                if (igTableSetColumnIndex(ic_q_qpt)) then
                   s = string(sys(isys)%c%vib%freq(i,win(w%idparent)%sc%iqpt_selected)*unitfactor,'f',&
                      length=9,decimal=digits,justify=ioj_right)
                   call iw_text(s)
                end if
             end do
          end if
          call igEndTable()
       end if ! igBeginTable (frequencies)
       call igEndGroup()

       ! suggested periodicity
       if (win(w%idparent)%sc%iqpt_selected > 0 .and..not.sys(isys)%c%ismolecule) then
          call igAlignTextToFramePadding()
          call iw_text("Suggested Periodicity:",highlight=.true.)
          do i = 1, 3
             if (abs(sys(isys)%c%vib%qpt(i,win(w%idparent)%sc%iqpt_selected)) > rational_approx_eps) then
                call rational_approx(sys(isys)%c%vib%qpt(i,win(w%idparent)%sc%iqpt_selected),q,r(i),rational_approx_eps)
             else
                r(i) = 1
             end if
          end do
          call iw_text("[" // string(r(1)) // " " // string(r(2)) // " " // string(r(3)) // "]",sameline=.true.)
          if (iw_button("Set",sameline=.true.)) then
             win(w%idparent)%sc%nc = int(r)
             win(w%idparent)%sc%forcebuildlists = .true.
          end if
          call iw_tooltip("Change the number of unit cells represented to the suggested value",ttshown)
       end if

       ! set initial value of animation to automatic
       if (win(w%idparent)%sc%iqpt_selected > 0 .and. win(w%idparent)%sc%ifreq_selected > 0 .and.&
          win(w%idparent)%sc%animation == 0 .and. .not.fset) then
          win(w%idparent)%sc%animation = 2
          win(w%idparent)%sc%anim_speed = anim_speed_default
          win(w%idparent)%sc%anim_amplitude = anim_amplitude_default
       end if

       ! animation radio buttons (no animation if no values selected)
       call iw_text("Animation",highlight=.true.)
       if (iw_radiobutton("None",int=win(w%idparent)%sc%animation,intval=0_c_int)) then
          win(w%idparent)%forcerender = .true.
       end if
       call iw_tooltip("Stop the animation",ttshown)
       if (iw_radiobutton("Automatic",int=win(w%idparent)%sc%animation,intval=2_c_int,sameline=.true.)) then
          win(w%idparent)%sc%anim_speed = anim_speed_default
          win(w%idparent)%sc%anim_amplitude = anim_amplitude_default
       end if
       call iw_tooltip("Animate the scene with atomic displacements corresponding to a periodic phase",ttshown)
       if (iw_radiobutton("Manual/Nudge Structure",int=win(w%idparent)%sc%animation,intval=1_c_int,sameline=.true.)) then
          win(w%idparent)%sc%anim_amplitude = 0._c_float
          win(w%idparent)%sc%anim_phase = 0._c_float
          win(w%idparent)%forcerender = .true.
       end if
       call iw_tooltip("Animate the scene using a manually set atomic displacement value",ttshown)

       if (win(w%idparent)%sc%animation == 1) then
          ! manual
          if (.not.sys(isys)%c%ismolecule) then
             ! crystals
             str1 = "Amplitude##amplitude" // c_null_char
             str2 = "%.2f" // c_null_char
             call igPushItemWidth(iw_calcwidth(5,1))
             if (igDragFloat(c_loc(str1),win(w%idparent)%sc%anim_amplitude,&
                0.01_c_float,0._c_float,anim_amplitude_max,c_loc(str2),ImGuiSliderFlags_AlwaysClamp))&
                win(w%idparent)%forcerender = .true.
             call igPopItemWidth()
             call iw_tooltip("Amplitude of the atomic displacements",ttshown)

             call igSameLine(0._c_float,-1._c_float)
             str1 = "Phase##phase" // c_null_char
             str2 = "%.3f" // c_null_char
             call igPushItemWidth(iw_calcwidth(6,1))
             if (igDragFloat(c_loc(str1),win(w%idparent)%sc%anim_phase,&
                0.001_c_float,-1._c_float,1._c_float,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)) &
                win(w%idparent)%forcerender = .true.
             call igPopItemWidth()
             call iw_tooltip("Phase for the atomic displacements along the chosen phonon normal mode",ttshown)
          else
             ! molecules
             str1 = "Displacement##amplitude" // c_null_char
             str2 = "%.2f" // c_null_char
             call igPushItemWidth(iw_calcwidth(5,1))
             if (igDragFloat(c_loc(str1),win(w%idparent)%sc%anim_amplitude,&
                0.01_c_float,-anim_amplitude_max,anim_amplitude_max,c_loc(str2),ImGuiSliderFlags_AlwaysClamp))&
                win(w%idparent)%forcerender = .true.
             call igPopItemWidth()
             call iw_tooltip("Extent of the atomic displacements",ttshown)
          end if

          ! create nudged system
          if (iw_button("Create Nudged System")) then
             if (allocated(seed)) deallocate(seed)
             allocate(seed(1))
             call sys(isys)%c%makeseed_nudged(seed(1),sys(isys)%c%vib%qpt(:,win(w%idparent)%sc%iqpt_selected),&
                sys(isys)%c%vib%vec(:,:,win(w%idparent)%sc%ifreq_selected,win(w%idparent)%sc%iqpt_selected),&
                real(win(w%idparent)%sc%anim_amplitude,8),real(win(w%idparent)%sc%anim_phase,8))
             call add_systems_from_seeds(1,seed)
             call launch_initialization_thread()
          end if
          call iw_tooltip("Create a new system with displaced atomic positions as shown in the view",ttshown)

       elseif (win(w%idparent)%sc%animation == 2) then
          ! automatic
          str1 = "Amplitude##amplitude" // c_null_char
          str2 = "%.2f" // c_null_char
          call igPushItemWidth(iw_calcwidth(5,1))
          ldum = igDragFloat(c_loc(str1),win(w%idparent)%sc%anim_amplitude,&
             0.01_c_float,0._c_float,anim_amplitude_max,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
          call igPopItemWidth()
          call iw_tooltip("Amplitude of the atomic displacements",ttshown)

          call igSameLine(0._c_float,-1._c_float)
          str1 = "Speed##speed" // c_null_char
          str2 = "%.2f" // c_null_char
          call igPushItemWidth(iw_calcwidth(5,1))
          if (igDragFloat(c_loc(str1),win(w%idparent)%sc%anim_speed,&
             0.02_c_float,0.0_c_float,anim_speed_max,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)) &
             win(w%idparent)%sc%timerefanimation = time
          call igPopItemWidth()
          call iw_tooltip("Speed of the atomic displacements",ttshown)
       end if
    end if ! vib_ok

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(5,1,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! final buttons: close
    if (iw_button("Close")) doquit = .true.
    call iw_tooltip("Close this window",ttshown)

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit = close the window
    if (doquit) then
       call w%end()
    end if

  end subroutine draw_vibrations

  !xx! private procedures

  !> Draw the atom selection table and return whether any item has
  !> been changed. sty = atom draw style, isys = current system.
  function atom_selection_widget(c,r) result(changed)
    use gui_main, only: g
    use scenes, only: draw_style_atom, draw_style_molecule
    use utils, only: iw_text, iw_combo_simple, iw_tooltip, iw_calcheight, iw_checkbox,&
       iw_clamp_color3, iw_calcwidth, iw_button, iw_coloredit3
    use crystalmod, only: crystal
    use global, only: iunit_ang, dunit0
    use tools_io, only: string, ioj_right
    type(crystal), intent(in) :: c
    type(representation), intent(inout) :: r
    logical :: changed

    logical :: domol, ldum, oksel
    logical(c_bool) :: ch
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: s, str1, str2, str3
    real*8 :: x0(3)
    type(ImVec2) :: sz0, szero
    integer :: ispc, i, iz, ncol, ic_next, isys
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    real(c_float) :: pos

    ! column IDs in atom table
    integer, parameter :: ic_id = 0
    integer, parameter :: ic_name = 1
    integer, parameter :: ic_z = 2
    integer, parameter :: ic_shown = 3
    integer, parameter :: ic_color = 4
    integer, parameter :: ic_radius = 5

    ! column IDs in molecule table
    integer, parameter :: im_id = 0
    integer, parameter :: im_nat = 1
    integer, parameter :: im_shown = 2
    integer, parameter :: im_color = 3
    integer, parameter :: im_radius = 4
    integer, parameter :: im_rest = 5

    type(ImVec4), parameter :: rgbsel = ImVec4(1._c_float,0.8_c_float,0.1_c_float,0.5_c_float)

    logical, save :: ttshown = .false. ! tooltip flag

    szero%x = 0
    szero%y = 0
    isys = r%id
    call iw_text("Atom Selection",highlight=.true.)

    ! selector and reset
    changed = .false.
    ch = .false.
    if (.not.c%ismolecule) then
       call iw_combo_simple("Atom types##atomtypeselection","Species"//c_null_char//&
          "Symmetry-unique" //c_null_char//"Cell"//c_null_char//c_null_char,&
          r%atom_style%type,changed=ch)
    else
       call iw_combo_simple("Atom types##atomtypeselection","Species"//c_null_char//"Atoms"//c_null_char//&
          c_null_char,r%atom_style%type,changed=ch)
    end if
    call iw_tooltip("Group atoms by these categories",ttshown)
    if (ch) then
       call r%reset_atom_style()
       changed = .true.
    end if

    ! whether to do the molecule column
    ncol = 7
    domol = (r%atom_style%type == 2 .or. (r%atom_style%type == 1 .and. c%ismolecule))
    if (domol) ncol = 8

    ! atom style table, for atoms
    oksel = .false.
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tableatomstyles" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(5,r%atom_style%ntype)+1,0,.false.)
    if (igBeginTable(c_loc(str1),ncol,flags,sz0,0._c_float)) then
       ! header setup
       str2 = "Id" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_id)

       str2 = "Atom" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_name)

       str2 = "Z " // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_z)

       str2 = "Show" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_shown)

       str2 = "Col" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_color)

       str2 = "Radius" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_radius)

       ic_next = ic_radius + 1
       if (domol) then
          str2 = "Mol" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_next)
          ic_next = ic_next + 1
       end if

       if (c%ismolecule) then
          str2 = "Coordinates (Å)" // c_null_char
       else
          str2 = "Coordinates (fractional)" // c_null_char
       end if
       flags = ImGuiTableColumnFlags_WidthStretch
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_next)
       call igTableSetupScrollFreeze(0, 1) ! top row always visible

       ! draw the header
       call igTableHeadersRow()
       call igTableSetColumnWidthAutoAll(igGetCurrentTable())

       ! start the clipper
       clipper = ImGuiListClipper_ImGuiListClipper()
       call ImGuiListClipper_Begin(clipper,r%atom_style%ntype,-1._c_float)

       ! draw the rows
       do while(ImGuiListClipper_Step(clipper))
          call c_f_pointer(clipper,clipper_f)
          do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
             call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
             if (r%atom_style%type == 0) then
                ! species
                ispc = i
             elseif (r%atom_style%type == 1) then
                ! nneq
                ispc = c%at(i)%is
             elseif (r%atom_style%type == 2) then
                ! ncel
                ispc = c%atcel(i)%is
             end if
             iz = c%spc(ispc)%z

             ! id
             if (igTableSetColumnIndex(ic_id)) then
                call igAlignTextToFramePadding()
                call iw_text(string(i))

                ! the selectable
                call igPushStyleColor_Vec4(ImGuiCol_HeaderHovered,rgbsel)
                pos = igGetCursorPosX()
                flags = ImGuiSelectableFlags_SpanAllColumns
                flags = ior(flags,ImGuiSelectableFlags_AllowItemOverlap)
                str2 = "##selectableatomtable" // string(isys) // c_null_char
                call igSameLine(0._c_float,0._c_float)
                ldum = igSelectable_Bool(c_loc(str2),.false._c_bool,flags,szero)
                call igSetCursorPosX(pos)
                if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                   if (igIsMouseHoveringRect(g%LastItemData%NavRect%min,g%LastItemData%NavRect%max,.false._c_bool)) then
                      win(win(r%idwin)%idparent)%sc%nselection = 1
                      win(win(r%idwin)%idparent)%sc%selection_type = r%atom_style%type
                      win(win(r%idwin)%idparent)%sc%selection(1) = i
                      win(win(r%idwin)%idparent)%forcerender = .true.
                      oksel = .true.
                   end if
                end if
                call igPopStyleColor(1)
             end if

             ! name
             if (igTableSetColumnIndex(ic_name)) call iw_text(string(c%spc(ispc)%name))

             ! Z
             if (igTableSetColumnIndex(ic_z)) call iw_text(string(iz))

             ! shown
             if (igTableSetColumnIndex(ic_shown)) then
                changed = changed .or. iw_checkbox("##tableshown" // string(i) ,r%atom_style%shown(i))
                call iw_tooltip("Toggle display of the atom/bond/label associated to this atom",ttshown)
             end if

             ! color
             if (igTableSetColumnIndex(ic_color)) then
                ch = iw_coloredit3("##tablecolor",r%atom_style%rgb(:,i),nolabel=.true.)
                call iw_tooltip("Atom color",ttshown)
                if (ch) then
                   r%atom_style%rgb(:,i) = min(r%atom_style%rgb(:,i),1._c_float)
                   r%atom_style%rgb(:,i) = max(r%atom_style%rgb(:,i),0._c_float)
                   changed = .true.
                end if
             end if

             ! radius
             if (igTableSetColumnIndex(ic_radius)) then
                str2 = "##tableradius" // string(i) // c_null_char
                str3 = "%.3f" // c_null_char
                call igPushItemWidth(iw_calcwidth(5,1))
                ch = igDragFloat(c_loc(str2),r%atom_style%rad(i),0.01_c_float,0._c_float,5._c_float,c_loc(str3),&
                   ImGuiSliderFlags_AlwaysClamp)
                call iw_tooltip("Radius of the sphere representing the atom",ttshown)
                if (ch) then
                   r%atom_style%rad(i) = max(r%atom_style%rad(i),0._c_float)
                   changed = .true.
                end if
                call igPopItemWidth()
             end if

             ! molecule
             ic_next = ic_radius + 1
             if (domol) then
                ! i is a complete list index in this case
                if (igTableSetColumnIndex(ic_next)) call iw_text(string(c%idatcelmol(1,i)))
                ic_next = ic_next + 1
             end if

             ! rest of info
             if (igTableSetColumnIndex(ic_next)) then
                s = ""
                if (r%atom_style%type > 0) then
                   if (c%ismolecule) then
                      x0 = (c%atcel(i)%r+c%molx0) * dunit0(iunit_ang)
                   elseif (r%atom_style%type == 1) then
                      x0 = c%at(i)%x
                   else
                      x0 = c%atcel(i)%x
                   endif
                   s = string(x0(1),'f',8,4,ioj_right) //" "// string(x0(2),'f',8,4,ioj_right) //" "//&
                      string(x0(3),'f',8,4,ioj_right)
                end if
                call iw_text(s)
             end if
          end do ! clipper indices
       end do ! clipper step

       ! end the clipper and the table
       call ImGuiListClipper_End(clipper)
       call igEndTable()
    end if

    ! reset selection
    if (.not.oksel) then
       if (win(win(r%idwin)%idparent)%sc%nselection > 0) &
          win(win(r%idwin)%idparent)%forcerender = .true.
       win(win(r%idwin)%idparent)%sc%nselection = 0
    end if

    ! style buttons: show/hide
    if (iw_button("Show All##showallatoms")) then
       r%atom_style%shown = .true.
       changed = .true.
    end if
    call iw_tooltip("Show all atoms/bonds/labels in the system",ttshown)
    if (iw_button("Hide All##hideallatoms",sameline=.true.)) then
       r%atom_style%shown = .false.
       changed = .true.
    end if
    call iw_tooltip("Hide all atoms/bonds/labels in the system",ttshown)
    if (iw_button("Toggle Show/Hide##toggleallatoms",sameline=.true.)) then
       do i = 1, r%atom_style%ntype
          r%atom_style%shown(i) = .not.r%atom_style%shown(i)
       end do
       changed = .true.
    end if
    call iw_tooltip("Toggle the show/hide status for all atoms/bonds/labels",ttshown)

    ! molecule selection
    ! initialized and more than one molecule
    if (r%mol_style%isinit .and. r%mol_style%ntype > 1) then
       call iw_text("Molecule Selection",highlight=.true.)

       ! molecule style table, for molecules
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_Resizable)
       flags = ior(flags,ImGuiTableFlags_Reorderable)
       flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       str1="##tablemolstyles" // c_null_char
       sz0%x = 0
       sz0%y = iw_calcheight(min(3,r%mol_style%ntype)+1,0,.false.)
       if (igBeginTable(c_loc(str1),6,flags,sz0,0._c_float)) then
          ! header setup
          str2 = "Id" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,im_id)

          str2 = "nat" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,im_nat)

          str2 = "Show" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,im_shown)

          str2 = "Tint" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,im_color)

          str2 = "Scale" // c_null_char
          flags = ImGuiTableColumnFlags_None
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,im_radius)

          if (c%ismolecule) then
             str2 = "Center of mass (Å)" // c_null_char
          else
             str2 = "Center of mass (fractional)" // c_null_char
          end if
          flags = ImGuiTableColumnFlags_WidthStretch
          call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,im_rest)
          call igTableSetupScrollFreeze(0, 1) ! top row always visible

          ! draw the header
          call igTableHeadersRow()
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())

          ! start the clipper
          clipper = ImGuiListClipper_ImGuiListClipper()
          call ImGuiListClipper_Begin(clipper,r%mol_style%ntype,-1._c_float)

          ! draw the rows
          do while(ImGuiListClipper_Step(clipper))
             call c_f_pointer(clipper,clipper_f)
             do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd

                call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)

                ! id
                if (igTableSetColumnIndex(im_id)) then
                   call igAlignTextToFramePadding()
                   call iw_text(string(i))
                end if

                ! nat
                if (igTableSetColumnIndex(im_nat)) call iw_text(string(c%mol(i)%nat))

                ! shown
                if (igTableSetColumnIndex(im_shown)) then
                   changed = changed .or. iw_checkbox("##tablemolshown" // string(i) ,r%mol_style%shown(i))
                   call iw_tooltip("Toggle display of all atoms in this molecule",ttshown)
                end if

                ! color
                if (igTableSetColumnIndex(im_color)) then
                   ch = iw_coloredit3("##tablemolcolor" // string(i),r%mol_style%tint_rgb(:,i),nolabel=.true.)
                   call iw_tooltip("Molecule color tint",ttshown)
                   if (ch) then
                      r%mol_style%tint_rgb(:,i) = min(r%mol_style%tint_rgb(:,i),1._c_float)
                      r%mol_style%tint_rgb(:,i) = max(r%mol_style%tint_rgb(:,i),0._c_float)
                      changed = .true.
                   end if
                end if

                ! radius
                if (igTableSetColumnIndex(im_radius)) then
                   str2 = "##tablemolradius" // string(i) // c_null_char
                   str3 = "%.3f" // c_null_char
                   call igPushItemWidth(iw_calcwidth(5,1))
                   ch = igDragFloat(c_loc(str2),r%mol_style%scale_rad(i),0.005_c_float,0._c_float,&
                      5._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
                   call iw_tooltip("Scale factor for the atomic radii in this molecule",ttshown)
                   if (ch) then
                      r%mol_style%scale_rad(i) = max(r%mol_style%scale_rad(i),0._c_float)
                      changed = .true.
                   end if
                   call igPopItemWidth()
                end if

                ! rest of info
                if (igTableSetColumnIndex(im_rest)) then
                   x0 = c%mol(i)%cmass(.false.)
                   if (c%ismolecule) then
                      x0 = (x0+c%molx0) * dunit0(iunit_ang)
                   else
                      x0 = c%c2x(x0)
                   endif
                   s = string(x0(1),'f',8,4,ioj_right) //" "// string(x0(2),'f',8,4,ioj_right) //" "//&
                      string(x0(3),'f',8,4,ioj_right)
                   call iw_text(s)
                end if
             end do ! clipper indices
          end do ! clipper step

          ! end the clipper and the table
          call ImGuiListClipper_End(clipper)
          call igEndTable()
       end if
    end if

  end function atom_selection_widget

end submodule view
