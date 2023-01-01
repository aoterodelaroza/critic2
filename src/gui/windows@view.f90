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

contains

  !xx! view

  !> Draw the view.
  module subroutine draw_view(w)
    use interfaces_opengl3
    use interfaces_cimgui
    use keybindings, only: is_bind_event, BIND_VIEW_INC_NCELL, BIND_VIEW_DEC_NCELL,&
       BIND_VIEW_ALIGN_A_AXIS, BIND_VIEW_ALIGN_B_AXIS, BIND_VIEW_ALIGN_C_AXIS,&
       BIND_VIEW_ALIGN_X_AXIS, BIND_VIEW_ALIGN_Y_AXIS, BIND_VIEW_ALIGN_Z_AXIS,&
       BIND_NAV_ROTATE, BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET,&
       BIND_NAV_MEASURE, bindnames, get_bind_keyname
    use scenes, only: reptype_atoms, reptype_unitcell, style_phong
    use utils, only: iw_calcheight, iw_calcwidth, iw_clamp_color3, iw_combo_simple
    use global, only: dunit0, iunit_ang
    use gui_main, only: sysc, sys, sys_init, nsys, g, io, fontsize
    use utils, only: iw_text, iw_button, iw_tooltip, iw_combo_simple
    use tools_io, only: string, ioj_right
    use param, only: newline
    class(window), intent(inout), target :: w

    integer :: i, j, nrep, id, ipad, is, icel, ineq
    type(ImVec2) :: szavail, sz0, sz1, szero, pos
    type(ImVec4) :: tintcol, bgcol, col
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    character(len=:), allocatable, target :: msg
    logical(c_bool) :: is_selected
    logical :: hover, chbuild, chrender, goodsys, ldum, ok, ismol
    logical(c_bool) :: isatom, isbond, islabels, isuc, ch
    integer(c_int) :: amax, flags, nc(3), ires, idx(4), viewtype
    real(c_float) :: scal, width, sqw, ratio, depth, rgba(4)
    real*8 :: x0(3)

    logical, save :: ttshown = .false. ! tooltip flag

    ! coordinate this with representation_menu in scenes module
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

    ! flags for shortcuts
    isatom = .false.
    isbond = .false.
    islabels = .false.
    isuc = .false.
    do i = 1, sysc(w%view_selected)%sc%nrep
       if (sysc(w%view_selected)%sc%rep(i)%isinit) then
          if (sysc(w%view_selected)%sc%rep(i)%type == reptype_atoms) then
             isatom = isatom .or. sysc(w%view_selected)%sc%rep(i)%atoms_display
             isbond = isbond .or. sysc(w%view_selected)%sc%rep(i)%bonds_display
             islabels = islabels .or. sysc(w%view_selected)%sc%rep(i)%labels_display
          elseif (sysc(w%view_selected)%sc%rep(i)%type == reptype_unitcell) then
             isuc = isuc .or. sysc(w%view_selected)%sc%rep(i)%shown
          end if
       end if
    end do

    ! whether the selected view system is a good system
    goodsys = (w%view_selected >= 1 .and. w%view_selected <= nsys)
    if (goodsys) goodsys = sysc(w%view_selected)%status == sys_init

    ! scene menu
    str1="##viewscenebutton" // c_null_char
    if (iw_button("Scene")) then
       call igOpenPopup_Str(c_loc(str1),ImGuiPopupFlags_None)
    end if
    if (goodsys) then
       if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_None)) then
          ! display shortcuts
          call iw_text("Display Shortcuts",highlight=.true.)
          str2 = "Atoms##atomsshortcut" // c_null_char
          if (igCheckbox(c_loc(str2),isatom)) then
             do i = 1, sysc(w%view_selected)%sc%nrep
                if (sysc(w%view_selected)%sc%rep(i)%isinit) then
                   if (sysc(w%view_selected)%sc%rep(i)%type == reptype_atoms) then
                      sysc(w%view_selected)%sc%rep(i)%atoms_display = isatom
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display atoms in all representations",ttshown)

          str2 = "Bonds##bondsshortcut" // c_null_char
          call igSameLine(0._c_float,-1._c_float)
          if (igCheckbox(c_loc(str2),isbond)) then
             do i = 1, sysc(w%view_selected)%sc%nrep
                if (sysc(w%view_selected)%sc%rep(i)%isinit) then
                   if (sysc(w%view_selected)%sc%rep(i)%type == reptype_atoms) then
                      sysc(w%view_selected)%sc%rep(i)%bonds_display = isbond
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display bonds in all representations",ttshown)

          str2 = "Labels##labelshortcut" // c_null_char
          call igSameLine(0._c_float,-1._c_float)
          if (igCheckbox(c_loc(str2),islabels)) then
             do i = 1, sysc(w%view_selected)%sc%nrep
                if (sysc(w%view_selected)%sc%rep(i)%isinit) then
                   if (sysc(w%view_selected)%sc%rep(i)%type == reptype_atoms) then
                      sysc(w%view_selected)%sc%rep(i)%labels_display = islabels
                   end if
                end if
             end do
             chbuild = .true.
          end if
          call iw_tooltip("Toggle display labels in all representations",ttshown)

          if (.not.sys(w%view_selected)%c%ismolecule) then
             str2 = "Unit Cell##ucshortcut" // c_null_char
             call igSameLine(0._c_float,-1._c_float)
             if (igCheckbox(c_loc(str2),isuc)) then
                do i = 1, sysc(w%view_selected)%sc%nrep
                   if (sysc(w%view_selected)%sc%rep(i)%isinit) then
                      if (sysc(w%view_selected)%sc%rep(i)%type == reptype_unitcell) then
                         sysc(w%view_selected)%sc%rep(i)%shown = isuc
                      end if
                   end if
                end do
                chbuild = .true.
             end if
             call iw_tooltip("Toggle display of unit cell representations",ttshown)
          end if

          ! number of cells selector
          if (.not.sys(w%view_selected)%c%ismolecule) then
             call igAlignTextToFramePadding()
             call iw_text("Periodicity",highlight=.true.)
             if (iw_button("Reset##periodicity",sameline=.true.)) then
                sysc(w%view_selected)%sc%nc = 1
                chbuild = .true.
             end if
             call iw_tooltip("Reset the number of cells to one",ttshown)

             ! calculate widths
             ipad = ceiling(log10(max(maxval(sysc(w%view_selected)%sc%nc),1) + 0.1))
             sqw = max(iw_calcwidth(1,1),igGetTextLineHeightWithSpacing())
             call igPushItemWidth(sqw)

             nc = sysc(w%view_selected)%sc%nc
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
             if (any(nc /= sysc(w%view_selected)%sc%nc)) then
                sysc(w%view_selected)%sc%nc = nc
                chbuild = .true.
             end if

             call igPopItemWidth()
          end if

          ! align view axis
          call iw_text("Axes Alignment",highlight=.true.)
          if (.not.sys(w%view_selected)%c%ismolecule) then
             if (iw_button("a")) then
                call sysc(w%view_selected)%sc%align_view_axis(1)
                chrender = .true.
             end if
             call iw_tooltip("Align the camera along the crystallographic a axis",ttshown)
             if (iw_button("b",sameline=.true.)) then
                call sysc(w%view_selected)%sc%align_view_axis(2)
                chrender = .true.
             end if
             call iw_tooltip("Align the camera along the crystallographic b axis",ttshown)
             if (iw_button("c",sameline=.true.)) then
                call sysc(w%view_selected)%sc%align_view_axis(3)
                chrender = .true.
             end if
             call iw_tooltip("Align the camera along the crystallographic c axis",ttshown)
          end if
          if (iw_button("x",sameline=.not.sys(w%view_selected)%c%ismolecule)) then
             call sysc(w%view_selected)%sc%align_view_axis(-1)
             chrender = .true.
          end if
          call iw_tooltip("Align the camera along the Cartesian x axis",ttshown)
          if (iw_button("y",sameline=.true.)) then
             call sysc(w%view_selected)%sc%align_view_axis(-2)
             chrender = .true.
          end if
          call iw_tooltip("Align the camera along the Cartesian y axis",ttshown)
          if (iw_button("z",sameline=.true.)) then
             call sysc(w%view_selected)%sc%align_view_axis(-3)
             chrender = .true.
          end if
          call iw_tooltip("Align the camera along the Cartesian z axis",ttshown)
          call igSameLine(0._c_float,-1._c_float)
          str2 = "Reset Distance##resetdistance" // c_null_char
          str3 = "%.2f" // c_null_char
          call igPushItemWidth(iw_calcwidth(5,1))
          ch = igDragFloat(c_loc(str2),sysc(w%view_selected)%sc%camresetdist,&
             0.01_c_float,0.1_c_float,5.0_c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
          if (ch) then
             call sysc(w%view_selected)%sc%reset()
             chrender = .true.
          end if
          call igPopItemWidth()
          call iw_tooltip("Ratio controlling distance from object when resetting camera",ttshown)

          ! object resolution
          call iw_text("Object Resolution",highlight=.true.)
          ires = sysc(w%view_selected)%sc%atom_res - 1
          call iw_combo_simple("Atoms##atomresselect","1: Carnby"//c_null_char//"2: Rough"//c_null_char//&
             "3: Normal"//c_null_char//"4: Good"//c_null_char//"5: Amazing"//c_nulL_char,ires)
          call iw_tooltip("Set the resolution of the spheres representing the atoms",ttshown)
          if (ires + 1 /= sysc(w%view_selected)%sc%atom_res) then
             sysc(w%view_selected)%sc%atom_res = ires + 1
             chrender = .true.
          end if
          ires = sysc(w%view_selected)%sc%bond_res - 1
          call iw_combo_simple("Bonds##bondresselect","1: Rough" // c_null_char // "2: Normal" // c_null_char //&
             "3: Good" // c_null_char,ires,sameline=.true.)
          call iw_tooltip("Set the resolution of the cylinders representing the bonds",ttshown)
          if (ires + 1 /= sysc(w%view_selected)%sc%bond_res) then
             sysc(w%view_selected)%sc%bond_res = ires + 1
             chrender = .true.
          end if

          ! scene style
          call iw_text("Appearance",highlight=.true.)
          call iw_combo_simple("Style##scenestyle","Simple"//c_null_char//"Realistic"&
             //c_null_char//c_null_char,sysc(w%view_selected)%sc%style,changed=ch)
          if (ch) then
             call sysc(w%view_selected)%sc%set_style_defaults()
             sysc(w%view_selected)%sc%forcebuildlists = .true.
          end if

          if (sysc(w%view_selected)%sc%style == style_phong) then
             !! phong-specific options !!
             call igPushItemWidth(iw_calcwidth(15,3))
             str2 = "Light Position" // c_null_char
             str3 = "%.1f" // c_null_char
             chrender = chrender .or. igDragFloat3(c_loc(str2),sysc(w%view_selected)%sc%lightpos,&
                0.5_c_float,-FLT_MAX,FLT_MAX,c_loc(str3),ImGuiSliderFlags_None)
             call iw_tooltip("Change the position of the light",ttshown)
             call igPopItemWidth()

             call igPushItemWidth(iw_calcwidth(5,1))
             str2 = "Ambient " // c_null_char
             str3 = "%.3f" // c_null_char
             chrender = chrender .or. igDragFloat(c_loc(str2),sysc(w%view_selected)%sc%ambient,&
                0.002_c_float,0._c_float,1._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the ambient light intensity",ttshown)
             call igSameLine(0._c_float,-1._c_float)
             str2 = "Diffuse" // c_null_char
             str3 = "%.3f" // c_null_char
             chrender = chrender .or. igDragFloat(c_loc(str2),sysc(w%view_selected)%sc%diffuse,&
                0.002_c_float,0._c_float,1._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the diffuse light intensity",ttshown)
             str2 = "Specular" // c_null_char
             str3 = "%.3f" // c_null_char
             chrender = chrender .or. igDragFloat(c_loc(str2),sysc(w%view_selected)%sc%specular,&
                0.002_c_float,0._c_float,1._c_float,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the specular light intensity",ttshown)
             call igSameLine(0._c_float,-1._c_float)
             str2 = "Shininess" // c_null_char
             str3 = "%.0f" // c_null_char
             chrender = chrender .or. igDragInt(c_loc(str2),sysc(w%view_selected)%sc%shininess,&
                1._c_float,0_c_int,256_c_int,c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Change the shininess of the light",ttshown)
             call igPopItemWidth()

             str2 = "Light" // c_null_char
             chrender = chrender .or. igColorEdit3(c_loc(str2),sysc(w%view_selected)%sc%lightcolor,&
                ImGuiColorEditFlags_NoInputs)
             call iw_tooltip("Change the color of the light",ttshown)
             call iw_clamp_color3(sysc(w%view_selected)%sc%lightcolor)
             call igSameLine(0._c_float,-1._c_float)
          end if

          ! background color
          str2 = "Background" // c_null_char
          chrender = chrender .or. igColorEdit3(c_loc(str2),sysc(w%view_selected)%sc%bgcolor,&
             ImGuiColorEditFlags_NoInputs)
          call iw_tooltip("Change the scene background color",ttshown)
          call iw_clamp_color3(sysc(w%view_selected)%sc%bgcolor)

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
                      sysc(i)%sc%nc = sysc(w%view_selected)%sc%nc
                   sysc(i)%sc%atom_res = sysc(w%view_selected)%sc%atom_res
                   sysc(i)%sc%bond_res = sysc(w%view_selected)%sc%bond_res
                   sysc(i)%sc%lightpos = sysc(w%view_selected)%sc%lightpos
                   sysc(i)%sc%ambient = sysc(w%view_selected)%sc%ambient
                   sysc(i)%sc%diffuse = sysc(w%view_selected)%sc%diffuse
                   sysc(i)%sc%specular = sysc(w%view_selected)%sc%specular
                   sysc(i)%sc%shininess = sysc(w%view_selected)%sc%shininess
                   sysc(i)%sc%lightcolor = sysc(w%view_selected)%sc%lightcolor
                   sysc(i)%sc%bgcolor = sysc(w%view_selected)%sc%bgcolor
                   sysc(i)%sc%camresetdist = sysc(w%view_selected)%sc%camresetdist
                end if
                call sysc(i)%sc%build_lists()
             end do
          end if
          call iw_tooltip("Apply these settings to all loaded systems",ttshown)
          if (iw_button("Reset",sameline=.true.,danger=.true.)) then
             call sysc(w%view_selected)%sc%init(w%view_selected)
             chbuild = .true.
          end if
          call iw_tooltip("Reset to the default settings",ttshown)

          call igEndPopup()
       end if
    end if
    call iw_tooltip("Change the view options",ttshown)

    ! gear menu
    str1="##viewgear" // c_null_char
    if (iw_button("Reps",sameline=.true.)) then
       call igOpenPopup_Str(c_loc(str1),ImGuiPopupFlags_None)
    end if
    if (goodsys) then
       if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_None)) then
          call igAlignTextToFramePadding()
          ! representations table
          call iw_text("Representations",highlight=.true.)

          ! add button
          ldum = iw_button("Add",sameline=.true.,popupcontext=ok,popupflags=ImGuiPopupFlags_MouseButtonLeft)
          if (ok) then
             str2 = "Atoms" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                id = sysc(w%view_selected)%sc%get_new_representation_id()
                call sysc(w%view_selected)%sc%rep(id)%init(w%view_selected,id,&
                   reptype_atoms,sysc(w%view_selected)%sc%style)
                chbuild = .true.
             end if
             call iw_tooltip("Represent atoms, bonds, and labels in the scene",ttshown)

             if (.not.sys(w%view_selected)%c%ismolecule) then
                str2 = "Unit Cell" // c_null_char
                if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                   id = sysc(w%view_selected)%sc%get_new_representation_id()
                   call sysc(w%view_selected)%sc%rep(id)%init(w%view_selected,id,&
                      reptype_unitcell,sysc(w%view_selected)%sc%style)
                   chbuild = .true.
                end if
                call iw_tooltip("Represent the unit cell",ttshown)
             end if
             call igEndPopup()
          end if
          call iw_tooltip("Add a representation to the view",ttshown)

          ! rest of the table
          str2 = "Representations##0,0" // c_null_char
          flags = ImGuiTableFlags_NoSavedSettings
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_NoBordersInBody)
          sz0%x = 0
          nrep = count(sysc(w%view_selected)%sc%rep(1:sysc(w%view_selected)%sc%nrep)%isinit)
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

             if (sysc(w%view_selected)%sc%representation_menu(w%id)) chbuild = .true.

             call igEndTable()
          end if

          call igEndPopup()
       end if
    end if
    call iw_tooltip("Add, remove, and modify representations",ttshown)

    ! update the draw lists and render
    if (chbuild) sysc(w%view_selected)%sc%forcebuildlists = .true.
    if (chrender .or. sysc(w%view_selected)%sc%forcebuildlists) w%forcerender = .true.

    ! save image
    if (iw_button("Export",sameline=.true.)) then
       if (w%idexportwin == 0) then
          w%idexportwin = stack_create_window(wintype_exportimage,.true.,idcaller=w%id)
       else
          call igSetWindowFocus_Str(c_loc(win(w%idexportwin)%name))
       end if
    end if
    call iw_tooltip("Export the current scene to an image file",ttshown)

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
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) then
                if (w%view_selected /= i) w%forcerender = .true.
                w%view_selected = i
                goodsys = .true.
             end if
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

    ! resize the render texture if not large enough
    amax = max(ceiling(max(szavail%x,szavail%y)),1)
    if (amax > w%FBOside) then
       amax = max(ceiling(1.5 * ceiling(max(szavail%x,szavail%y))),1)
       call w%delete_texture_view()
       call w%create_texture_view(amax)
       w%forcerender = .true.
    end if

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
    sysc(w%view_selected)%sc%camratio = min(ratio,2.5_c_float)

    ! render the image to the texture, if requested
    if (w%forcerender) then
       ! render to the draw framebuffer
       call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
       call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
       if (goodsys) then
          call glClearColor(sysc(w%view_selected)%sc%bgcolor(1),sysc(w%view_selected)%sc%bgcolor(2),&
             sysc(w%view_selected)%sc%bgcolor(3),1._c_float)
       else
          call glClearColor(0._c_float,0._c_float,0._c_float,0._c_float)
       end if
       call glClear(ior(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT),GL_STENCIL_BUFFER_BIT))
       if (goodsys) &
          call sysc(w%view_selected)%sc%render()
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! render to the pick frame buffer
       call glBindFramebuffer(GL_FRAMEBUFFER, w%FBOpick)
       call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
       call glClearColor(0._c_float,0._c_float,0._c_float,0._c_float)
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
       if (goodsys) call sysc(w%view_selected)%sc%renderpick()
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
    ldum = igImageButton(w%FBOtex, szavail, sz0, sz1, 0_c_int, bgcol, tintcol)
    call igPopStyleColor(3)

    ! get hover, image rectangle coordinates, and atom idx
    idx = 0
    hover = igIsItemHovered(ImGuiHoveredFlags_None) .and. goodsys
    call igGetItemRectMin(w%v_rmin)
    call igGetItemRectMax(w%v_rmax)
    if (hover) then
       call igGetMousePos(pos)
       call w%mousepos_to_texpos(pos)
       call w%getpixel(w%FBOpick,pos,depth,rgba)
       idx = transfer(rgba,idx)
    end if

    ! mode selection
    viewtype = 0
    call iw_combo_simple("##viewmode","Navigate" // c_null_char,viewtype)
    msg = trim(get_bind_keyname(BIND_NAV_ROTATE)) // ": " // trim(bindnames(BIND_NAV_ROTATE)) // newline
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
       ineq = sys(w%view_selected)%c%atcel(icel)%is
       ismol = sys(w%view_selected)%c%ismolecule

       msg = trim(sys(w%view_selected)%c%spc(is)%name)
       if (.not.ismol) then
          x0 = (sys(w%view_selected)%c%atcel(icel)%r+sys(w%view_selected)%c%molx0) * dunit0(iunit_ang)

          msg = trim(msg) // " [cellid=" // string(icel) // "+(" // string(idx(2)) // "," // string(idx(3)) //&
             "," // string(idx(4)) // "),nneqid=" // string(ineq) // ",wyckoff=" // &
             string(sys(w%view_selected)%c%at(ineq)%mult) // string(sys(w%view_selected)%c%at(ineq)%wyc)
          if (sys(w%view_selected)%c%nmol > 1) &
             msg = msg // ",molid=" // string(sys(w%view_selected)%c%idatcelmol(icel))
          msg = msg // "] " //&
             string(x0(1),'f',decimal=4) //" "// string(x0(2),'f',decimal=4) //" "//&
             string(x0(3),'f',decimal=4) // " (frac)"
       else
          x0 = sys(w%view_selected)%c%atcel(icel)%x

          msg = trim(msg) // " [id=" // string(icel)
          if (sys(w%view_selected)%c%nmol > 1) &
             msg = msg // ",molid=" // string(sys(w%view_selected)%c%idatcelmol(icel))
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
    if (goodsys) then
       if (.not.sys(w%view_selected)%c%ismolecule) then
          if (is_bind_event(BIND_VIEW_INC_NCELL)) then
             sysc(w%view_selected)%sc%nc = sysc(w%view_selected)%sc%nc + 1
             sysc(w%view_selected)%sc%forcebuildlists = .true.
          elseif (is_bind_event(BIND_VIEW_DEC_NCELL)) then
             sysc(w%view_selected)%sc%nc = max(sysc(w%view_selected)%sc%nc - 1,1)
             sysc(w%view_selected)%sc%forcebuildlists = .true.
          end if
          if (is_bind_event(BIND_VIEW_ALIGN_A_AXIS)) then
             call sysc(w%view_selected)%sc%align_view_axis(1)
             w%forcerender = .true.
          end if
          if (is_bind_event(BIND_VIEW_ALIGN_B_AXIS)) then
             call sysc(w%view_selected)%sc%align_view_axis(2)
             w%forcerender = .true.
          end if
          if (is_bind_event(BIND_VIEW_ALIGN_C_AXIS)) then
             call sysc(w%view_selected)%sc%align_view_axis(3)
             w%forcerender = .true.
          end if
       end if
       if (is_bind_event(BIND_VIEW_ALIGN_X_AXIS)) then
          call sysc(w%view_selected)%sc%align_view_axis(-1)
          w%forcerender = .true.
       end if
       if (is_bind_event(BIND_VIEW_ALIGN_Y_AXIS)) then
          call sysc(w%view_selected)%sc%align_view_axis(-2)
          w%forcerender = .true.
       end if
       if (is_bind_event(BIND_VIEW_ALIGN_Z_AXIS)) then
          call sysc(w%view_selected)%sc%align_view_axis(-3)
          w%forcerender = .true.
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
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
    call glBindTexture(GL_TEXTURE_2D, 0)

    call glBindTexture(GL_TEXTURE_2D, w%FBOrgba)
    call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, atex, atex, 0, GL_RGBA, GL_FLOAT, c_null_ptr)
    call glBindTexture(GL_TEXTURE_2D, 0)

    ! render buffers
    call glBindRenderbuffer(GL_RENDERBUFFER, w%FBOdepth)
    call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, atex, atex)
    call glBindRenderbuffer(GL_RENDERBUFFER, 0)
    call glBindRenderbuffer(GL_RENDERBUFFER, w%FBOdepthp)
    call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, atex, atex)
    call glBindRenderbuffer(GL_RENDERBUFFER, 0)

    ! frame buffers
    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
    call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, w%FBOtex, 0)
    call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, w%FBOdepth)
    call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_STENCIL_ATTACHMENT, GL_RENDERBUFFER, w%FBOdepth)
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
    call glClear(ior(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT),GL_STENCIL_BUFFER_BIT))
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
    use gui_main, only: nsys
    class(window), intent(inout), target :: w
    integer, intent(in) :: isys

    if (isys < 1 .or. isys > nsys) return
    if (w%view_selected == isys) return
    w%view_selected = isys
    w%forcerender = .true.

  end subroutine select_view

  !> Process the mouse events in the view window
  module subroutine process_events_view(w,hover,idx)
    use interfaces_cimgui
    use scenes, only: scene, min_zoom, max_zoom
    use utils, only: translate, rotate, mult, invmult
    use tools_math, only: cross_cfloat, matinv_cfloat
    use keybindings, only: is_bind_event, is_bind_mousescroll, BIND_NAV_ROTATE,&
       BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET, BIND_NAV_MEASURE,&
       BIND_CLOSE_FOCUSED_DIALOG
    use gui_main, only: io, nsys, sysc
    class(window), intent(inout), target :: w
    logical, intent(in) :: hover
    integer(c_int), intent(in) :: idx(4)

    type(ImVec2) :: texpos, mousepos
    real(c_float) :: ratio, pos3(3), vnew(3), vold(3), axis(3), lax
    real(c_float) :: mpos2(2), ang, xc(3)
    type(scene), pointer :: sc

    integer, parameter :: ilock_no = 1
    integer, parameter :: ilock_left = 2
    integer, parameter :: ilock_right = 3
    integer, parameter :: ilock_scroll = 4

    real(c_float), parameter :: mousesens_zoom0 = 0.15_c_float
    real(c_float), parameter :: mousesens_rot0 = 3._c_float

    type(ImVec2), save :: mposlast
    real(c_float), save :: mpos0_r(3), mpos0_l(3), cpos0_l(3)
    real(c_float), save :: oldview(4,4)
    real(c_float), save :: mpos0_s
    integer, save :: ilock = ilock_no

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! only process if there is an associated system is viewed and scene is initialized
    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    sc => sysc(w%view_selected)%sc
    if (.not.sc%isinit) return

    ! process mode-specific events
    if (w%view_mousebehavior == MB_Navigation) then
       call igGetMousePos(mousepos)
       texpos = mousepos

       ! transform to the texture pos
       call w%mousepos_to_texpos(texpos)

       ! Zoom. There are two behaviors: mouse scroll and hold key and
       ! translate mouse
       ratio = 0._c_float
       if (hover.and.(ilock == ilock_no .or. ilock == ilock_scroll).and. is_bind_event(BIND_NAV_ZOOM,.false.)) then
          if (is_bind_mousescroll(BIND_NAV_ZOOM)) then
             ! mouse scroll
             ratio = mousesens_zoom0 * io%MouseWheel
          else
             ! keys
             mpos0_s = mousepos%y
             ilock = ilock_scroll
          end if
       elseif (ilock == ilock_scroll) then
          if (is_bind_event(BIND_NAV_ZOOM,.true.)) then
             ! 10/a to make it adimensional
             ratio = mousesens_zoom0 * (mpos0_s-mousepos%y) * (10._c_float / w%FBOside)
             mpos0_s = mousepos%y
          else
             ilock = ilock_no
          end if
       end if
       if (ratio /= 0._c_float) then
          ratio = min(max(ratio,-0.99999_c_float),0.9999_c_float)

          xc = mult(sc%world,sc%scenecenter)
          pos3 = sc%campos - xc
          pos3 = pos3 - ratio * pos3
          if (norm2(pos3) < min_zoom) &
             pos3 = pos3 / norm2(pos3) * min_zoom
          if (norm2(pos3) > max_zoom * sc%scenerad) &
             pos3 = pos3 / norm2(pos3) * (max_zoom * sc%scenerad)
          sc%campos = xc + pos3

          call sc%update_view_matrix()
          call sc%update_projection_matrix()
          w%forcerender = .true.
       end if

       ! drag
       if (hover .and. is_bind_event(BIND_NAV_TRANSLATE,.false.) .and. (ilock == ilock_no .or. ilock == ilock_right)) then
          mpos0_r = (/texpos%x,texpos%y,1._c_float/)

          ! save the current view matrix
          oldview = sc%view

          ilock = ilock_right
          mposlast = mousepos
       elseif (ilock == ilock_right) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (is_bind_event(BIND_NAV_TRANSLATE,.true.)) then
             if (mousepos%x /= mposlast%x .or. mousepos%y /= mposlast%y) then
                vnew = (/texpos%x,texpos%y,mpos0_r(3)/)
                call w%texpos_to_view(vnew)
                vold = mpos0_r
                call w%texpos_to_view(vold)

                xc = vold - vnew
                xc = invmult(oldview,xc)
                sc%campos = xc
                call sc%update_view_matrix()
                w%forcerender = .true.
             end if
          else
             ilock = ilock_no
          end if
       end if

       ! rotate
       if (hover .and. is_bind_event(BIND_NAV_ROTATE,.false.) .and. (ilock == ilock_no .or. ilock == ilock_left)) then
          mpos0_l = (/texpos%x, texpos%y, 0._c_float/)
          cpos0_l = mpos0_l
          call w%texpos_to_view(cpos0_l)
          ilock = ilock_left
       elseif (ilock == ilock_left) then
          call igSetMouseCursor(ImGuiMouseCursor_Hand)
          if (is_bind_event(BIND_NAV_ROTATE,.true.)) then
             if (texpos%x /= mpos0_l(1) .or. texpos%y /= mpos0_l(2)) then
                vnew = (/texpos%x,texpos%y,mpos0_l(3)/)
                call w%texpos_to_view(vnew)
                pos3 = (/0._c_float,0._c_float,1._c_float/)
                axis = cross_cfloat(pos3,vnew - cpos0_l)
                lax = norm2(axis)
                if (lax > 1e-10_c_float) then
                   axis = axis / lax
                   axis = invmult(sc%world,axis,notrans=.true.)
                   mpos2(1) = texpos%x - mpos0_l(1)
                   mpos2(2) = texpos%y - mpos0_l(2)
                   ang = 2._c_float * norm2(mpos2) * mousesens_rot0 / w%FBOside
                   sc%world = translate(sc%world,sc%scenecenter)
                   sc%world = rotate(sc%world,ang,axis)
                   sc%world = translate(sc%world,-sc%scenecenter)

                   w%forcerender = .true.
                end if
                mpos0_l = (/texpos%x, texpos%y, 0._c_float/)
                cpos0_l = mpos0_l
                call w%texpos_to_view(cpos0_l)
             end if
          else
             ilock = ilock_no
          end if
       end if

       ! reset the view
       if (hover .and. is_bind_event(BIND_NAV_RESET,.false.)) then
          call sc%reset()
          w%forcerender = .true.
       end if

       ! atom selection
       if (hover .and. is_bind_event(BIND_NAV_MEASURE)) then
          call sc%select_atom(idx)
          w%forcerender = .true.
       end if
       if (hover .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) then
          call sc%select_atom((/0,0,0,0/))
          w%forcerender = .true.
       end if
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      mposlast%x = 0._c_float
      mposlast%y = 0._c_float
      mpos0_r = 0._c_float
      mpos0_l = 0._c_float
      oldview = 0._c_float
      cpos0_l = 0._c_float
      mpos0_s = 0._c_float
      ilock = ilock_no
    end subroutine init_state

  end subroutine process_events_view

  module subroutine draw_selection_tooltip(w,idx)
    use interfaces_cimgui
    use utils, only: iw_text
    use gui_main, only: sys, sysc, fontsize
    use tools_io, only: string
    use tools_math, only: cross
    use param, only: bohrtoa, pi
    class(window), intent(inout), target :: w
    integer(c_int), intent(in) :: idx(4)

    character(kind=c_char,len=:), allocatable, target :: msg

    integer :: nmsel, i
    integer :: msel(4,4)
    type(ImVec4) :: col
    integer :: idx1(4), idx2(4), idx3(4), idx4(4)
    real*8 :: x0(3), x1(3), x2(3), d, d1, d2, ang, n0(3), n1(3)


    real(c_float), parameter :: rgbsel(4,4) = reshape((/&
       1._c_float,  0.4_c_float, 0.4_c_float, 1._c_float,&
       0.4_c_float, 1._c_float,  0.4_c_float, 1._c_float,&
       0.4_c_float, 0.4_c_float, 1._c_float, 1._c_float,&
       0.9_c_float, 0.7_c_float, 0.4_c_float, 1._c_float/),shape(rgbsel))

    ! check if the tooltip is needed
    nmsel = sysc(w%view_selected)%sc%nmsel
    if (nmsel == 0) return
    msel = sysc(w%view_selected)%sc%msel
    if (nmsel == 1 .and. (idx(1) == 0 .or. all(idx == msel(:,1)))) return

    ! start tooltip and header
    call igBeginTooltip()
    call igPushTextWrapPos(60._c_float * fontsize%x)
    call iw_text("Distance (d), angle (α), dihedral (φ)")

    ! distance 1-2
    idx1 = msel(:,1)
    if (nmsel == 1) then
       idx2 = idx
       if (idx(1) == 0) goto 999
    else
       idx2 = msel(:,2)
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
       idx1 = msel(:,2)
       if (nmsel == 2) then
          idx2 = idx
          if (idx(1) == 0) goto 999
          if (all(idx == msel(:,2)) .or. all(idx == msel(:,1))) goto 999
       else
          idx2 = msel(:,3)
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
       idx3 = msel(:,1)
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
       idx1 = msel(:,3)
       if (nmsel == 3) then
          idx2 = idx
          if (idx(1) == 0) goto 999
          if (all(idx == msel(:,3)) .or. all(idx == msel(:,2)) .or. all(idx == msel(:,1))) goto 999
       else
          idx2 = msel(:,4)
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
       idx3 = msel(:,2)
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
       idx4 = msel(:,1)
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
  module subroutine getpixel(w,fb,pos,depth,rgba)
    use interfaces_opengl3
    class(window), intent(inout), target :: w
    integer(c_int), intent(in) :: fb
    type(ImVec2), intent(inout) :: pos
    real(c_float), intent(out), optional :: depth
    real(c_float), intent(out), optional :: rgba(4)

    real(c_float), target :: depth_, rgba_(4)

    call glBindFramebuffer(GL_FRAMEBUFFER, fb)
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
    use gui_main, only: sysc, nsys
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    type(scene), pointer :: sc

    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    sc => sysc(w%view_selected)%sc
    if (.not.sc%isinit) return

    pos = project(pos,eye4,sc%projection,w%FBOside)

  end subroutine view_to_texpos

  !> Transform texture position (x,y) plus depth (z) to view
  !> coordinates.
  module subroutine texpos_to_view(w,pos)
    use utils, only: unproject
    use scenes, only: scene
    use gui_main, only: sysc, nsys
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    type(scene), pointer :: sc

    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    sc => sysc(w%view_selected)%sc
    if (.not.sc%isinit) return

    pos = unproject(pos,eye4,sc%projection,w%FBOside)

  end subroutine texpos_to_view

  !> Transform from world coordinates to texture position (x,y)
  !> plus depth (z).
  module subroutine world_to_texpos(w,pos)
    use utils, only: project
    use scenes, only: scene
    use gui_main, only: sysc, nsys
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    type(scene), pointer :: sc

    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    sc => sysc(w%view_selected)%sc
    if (.not.sc%isinit) return

    pos = project(pos,sc%view,sc%projection,w%FBOside)

  end subroutine world_to_texpos

  !> Transform texture position (x,y) plus depth (z) to world
  !> coordinates.
  module subroutine texpos_to_world(w,pos)
    use utils, only: unproject
    use scenes, only: scene
    use gui_main, only: sysc, nsys
    class(window), intent(inout), target :: w
    real(c_float), intent(inout) :: pos(3)

    type(scene), pointer :: sc

    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    sc => sysc(w%view_selected)%sc
    if (.not.sc%isinit) return

    pos = unproject(pos,sc%view,sc%projection,w%FBOside)

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
       iw_radiobutton, iw_calcheight
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
    if (.not.doquit) doquit = win(w%idparent)%type /= wintype_view

    if (.not.doquit) then
       ! whether the rep has changed
       changed = .false.

       !!! name and type block
       call iw_text("Name and Type",highlight=.true.)

       ! the representation type
       itype = w%rep%type - 1
       call iw_combo_simple("##reptype","Atoms" // c_null_char // "Unit cell" // c_null_char,itype)
       if (w%rep%type /= itype + 1) changed = .true.
       w%rep%type = itype + 1
       call iw_tooltip("Type of representation",ttshown)

       ! name text input
       str1 = "##nametextinput"
       txtinp = trim(adjustl(w%rep%name)) // c_null_char
       call igSameLine(0._c_float,-1._c_float)
       call igPushItemWidth(iw_calcwidth(30,1))
       if (igInputText(c_loc(str1),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_None,c_null_ptr,c_null_ptr)) then
          ll = index(txtinp,c_null_char)
          w%rep%name = txtinp(1:ll-1)
       end if
       call igPopItemWidth()
       call iw_tooltip("Name of this representation",ttshown)

       ! shown checkbox
       call igSameLine(0._c_float,-1._c_float)
       str2 = "Show" // c_null_char
       changed = changed .or. igCheckbox(c_loc(str2),w%rep%shown)
       call iw_tooltip("Toggle show/hide this representation",ttshown)

       ! type-dependent items
       if (w%rep%type == reptype_atoms) then
          changed = changed .or. w%draw_editrep_atoms(ttshown)
       elseif (w%rep%type == reptype_unitcell) then
          changed = changed .or. w%draw_editrep_unitcell(ttshown)
       end if

       ! rebuild draw lists if necessary
       if (changed) sysc(win(w%idparent)%view_selected)%sc%forcebuildlists = .true.

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
          call w%rep%init(w%rep%id,w%rep%idrep,itype,sysc(win(w%idparent)%view_selected)%sc%style)
          w%rep%name = str2
          w%rep%shown = lshown
          sysc(win(w%idparent)%view_selected)%sc%forcebuildlists = .true.
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

  !> Draw the editrep window, atoms class. Returns true if the scene needs
  !> rendering again. ttshown = the tooltip flag.
  module function draw_editrep_atoms(w,ttshown) result(changed)
    use global, only: dunit0, iunit_ang
    use scenes, only: representation
    use gui_main, only: sys, g, ColorHighlightText
    use utils, only: iw_text, iw_tooltip, iw_combo_simple, iw_button, iw_calcwidth,&
       iw_radiobutton, iw_calcheight, iw_clamp_color3
    use tools_io, only: string, ioj_right, ioj_left
    use param, only: atmcov, atmvdw, jmlcol, jmlcol2, newline
    class(window), intent(inout), target :: w
    logical, intent(inout) :: ttshown
    logical(c_bool) :: changed

    integer :: ispc, isys, iz, ll, ipad
    character(kind=c_char,len=1024), target :: txtinp
    character(len=:), allocatable :: s
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    type(ImVec2) :: sz0
    real*8 :: x0(3), res
    logical(c_bool) :: ch, ldum
    integer(c_int) :: flags, nc(3), lst
    real(c_float) :: sqw
    integer :: i

    integer(c_int), parameter :: lsttrans(0:6) = (/0,1,1,1,2,3,4/)
    integer(c_int), parameter :: lsttransi(0:4) = (/0,1,4,5,6/)

    integer, parameter :: ic_id = 0
    integer, parameter :: ic_name = 1
    integer, parameter :: ic_z = 2
    integer, parameter :: ic_shown = 3
    integer, parameter :: ic_color = 4
    integer, parameter :: ic_radius = 5
    integer, parameter :: ic_rest = 6

    ! initialize
    changed = .false.
    isys = w%isys

    ! filter
    call igAlignTextToFramePadding()
    call iw_text("Filter",highlight=.true.)
    call iw_text("(?)",sameline=.true.)
    call iw_tooltip("Show the atom if the filter expression evaluates to non-zero (true) at the atomic position. "&
       &"Structural variables are very useful for filters. Examples:"//newline//&
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
    str1 = "##filtertext"
    txtinp = trim(adjustl(w%rep%filter)) // c_null_char
    if (igInputText(c_loc(str1),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_EnterReturnsTrue,&
       c_null_ptr,c_null_ptr)) then
       ll = index(txtinp,c_null_char)
       w%rep%filter = txtinp(1:ll-1)

       ! test the filter
       if (sys(isys)%c%ncel > 0) then
          x0 = sys(isys)%c%atcel(1)%r
       else
          x0 = 0d0
       end if
       res = sys(isys)%eval(w%rep%filter,.false.,w%rep%goodfilter,x0)

       if (w%rep%goodfilter) changed = .true.
    end if
    if (len_trim(w%rep%filter) == 0) w%rep%goodfilter = .true.
    call iw_tooltip("Apply this filter to the atoms in the system. Atoms are represented if non-zero.",&
       ttshown)
    if (iw_button("Clear",sameline=.true.)) then
       w%rep%filter = ""
       w%rep%goodfilter = .true.
       changed = .true.
    end if
    call iw_tooltip("Clear the filter",ttshown)
    if (.not.w%rep%goodfilter) &
       call iw_text("Error in the filter expression.",danger=.true.)

    ! periodicity
    if (.not.sys(isys)%c%ismolecule) then
       call igAlignTextToFramePadding()
       call iw_text("Periodicity",highlight=.true.)

       ! radio buttons for the periodicity type
       changed = changed .or. iw_radiobutton("None",int=w%rep%pertype,intval=0_c_int,sameline=.true.)
       call iw_tooltip("Cell not repeated for this representation",ttshown)
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

       ! checkbox for molecular motif
       if (sys(isys)%c%ismol3d) then
          str2 = "Show connected molecules" // c_null_char
          changed = changed .or. igCheckbox(c_loc(str2),w%rep%onemotif)
          call iw_tooltip("Translate atoms to display whole molecules",ttshown)
       end if

       ! checkbox for border
       str2 = "Show atoms at cell edges" // c_null_char
       changed = changed .or. igCheckbox(c_loc(str2),w%rep%border)
       call iw_tooltip("Display atoms near the unit cell edges",ttshown)
    end if

    ! global display control
    call iw_text("Display and Atom Selection",highlight=.true.)

    ! origin of the atoms
    call igPushItemWidth(iw_calcwidth(21,3))
    if (sys(isys)%c%ismolecule) then
       str1 = "Origin shift (Å)##originatom" // c_null_char
       str2 = "%.4f" // c_null_char
       changed = changed .or. igDragFloat3(c_loc(str1),w%rep%origin,&
          0.01_c_float,-FLT_MAX,FLT_MAX,c_loc(str2),ImGuiSliderFlags_None)
    else
       str1 = "Origin shift (fractional coord.)##originatom" // c_null_char
       str2 = "%.5f" // c_null_char
       changed = changed .or. igDragFloat3(c_loc(str1),w%rep%origin,&
          0.001_c_float,-FLT_MAX,FLT_MAX,c_loc(str2),ImGuiSliderFlags_None)
    end if
    call iw_tooltip("Coordinates for the origin shift of the atoms/bonds/labels",ttshown)
    call igPopItemWidth()

    ! selector and reset
    ch = .false.
    if (.not.sys(isys)%c%ismolecule) then
       call iw_combo_simple("Atom types##atomtypeselection","Species"//c_null_char//&
          "Symmetry-unique" //c_null_char//"Cell"//c_null_char//c_null_char,&
          w%rep%atom_style_type,changed=ch)
    else
       call iw_combo_simple("Atom types##atomtypeselection","Species"//c_null_char//"Atoms"//c_null_char//&
          c_null_char,w%rep%atom_style_type,changed=ch)
    end if
    call iw_tooltip("Group atoms in the representation by these categories",ttshown)
    if (ch) then
       call w%rep%reset_atom_style()
       changed = .true.
    end if

    ! atom style table, for atoms
    flags = ImGuiTableFlags_None
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
    flags = ior(flags,ImGuiTableFlags_Borders)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    str1="##tableatomstyles" // c_null_char
    sz0%x = 0
    sz0%y = iw_calcheight(min(5,w%rep%natom_style)+1,0,.false.)
    if (igBeginTable(c_loc(str1),8,flags,sz0,0._c_float)) then
       ! header setup
       str2 = "Id" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_id)

       str2 = "Atom" // c_null_char
       flags = ImGuiTableColumnFlags_None
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_name)

       str2 = "Z" // c_null_char
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

       if (sys(isys)%c%ismolecule) then
          str2 = "Coordinates (Å)" // c_null_char
       else
          str2 = "Coordinates (fractional)" // c_null_char
       end if
       flags = ImGuiTableColumnFlags_WidthStretch
       call igTableSetupColumn(c_loc(str2),flags,0.0_c_float,ic_rest)
       call igTableSetupScrollFreeze(0, 1) ! top row always visible

       ! draw the header
       call igTableHeadersRow()
       call igTableSetColumnWidthAutoAll(igGetCurrentTable())

       ! draw the rows
       do i = 1, w%rep%natom_style
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
          if (w%rep%atom_style_type == 0) then
             ! species
             ispc = i
          elseif (w%rep%atom_style_type == 1) then
             ! nneq
             ispc = sys(isys)%c%at(i)%is
          elseif (w%rep%atom_style_type == 2) then
             ! ncel
             ispc = sys(isys)%c%atcel(i)%is
          end if
          iz = sys(isys)%c%spc(ispc)%z

          ! id
          if (igTableSetColumnIndex(ic_id)) call iw_text(string(i))

          ! name
          if (igTableSetColumnIndex(ic_name)) call iw_text(string(sys(isys)%c%spc(ispc)%name))

          ! Z
          if (igTableSetColumnIndex(ic_z)) call iw_text(string(iz))

          ! shown
          if (igTableSetColumnIndex(ic_shown)) then
             str2 = "##tableshown" // string(i) // c_null_char
             changed = changed .or. igCheckbox(c_loc(str2),w%rep%atom_style(i)%shown)
             call iw_tooltip("Toggle display of the atom/bond/label associated to this atom",ttshown)
          end if

          ! color
          if (igTableSetColumnIndex(ic_color)) then
             str2 = "##tablecolor" // string(i) // c_null_char
             flags = ior(ImGuiColorEditFlags_NoInputs,ImGuiColorEditFlags_NoLabel)
             ch = igColorEdit3(c_loc(str2),w%rep%atom_style(i)%rgb,flags)
             call iw_tooltip("Atom color",ttshown)
             call iw_clamp_color3(w%rep%atom_style(i)%rgb)
             if (ch) then
                w%rep%atom_style(i)%rgb = min(w%rep%atom_style(i)%rgb,1._c_float)
                w%rep%atom_style(i)%rgb = max(w%rep%atom_style(i)%rgb,0._c_float)
                changed = .true.
             end if
          end if

          ! radius
          if (igTableSetColumnIndex(ic_radius)) then
             str2 = "##tableradius" // string(i) // c_null_char
             str3 = "%.3f" // c_null_char
             call igPushItemWidth(iw_calcwidth(5,1))
             ch = igDragFloat(c_loc(str2),w%rep%atom_style(i)%rad,0.01_c_float,0._c_float,5._c_float,c_loc(str3),&
                ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Radius of the sphere representing the atom",ttshown)
             if (ch) then
                w%rep%atom_style(i)%rad = max(w%rep%atom_style(i)%rad,0._c_float)
                changed = .true.
             end if
             call igPopItemWidth()
          end if

          ! rest of info
          if (igTableSetColumnIndex(ic_rest)) then
             s = ""
             if (w%rep%atom_style_type > 0) then
                if (sys(isys)%c%ismolecule) then
                   x0 = (sys(isys)%c%atcel(i)%r+sys(isys)%c%molx0) * dunit0(iunit_ang)
                elseif (w%rep%atom_style_type == 1) then
                   x0 = sys(isys)%c%at(i)%x
                elseif (w%rep%atom_style_type == 2) then
                   x0 = sys(isys)%c%atcel(i)%x
                endif
                if (w%rep%atom_style_type == 2 .or. (w%rep%atom_style_type == 1 .and. sys(isys)%c%ismolecule)) then
                   s = "mol="//string(sys(isys)%c%idatcelmol(i),2,ioj_left)
                end if
                s = s // string(x0(1),'f',8,4,ioj_right) //" "// string(x0(2),'f',8,4,ioj_right) //" "//&
                   string(x0(3),'f',8,4,ioj_right)
             end if
             call iw_text(s)
          end if
       end do
       call igEndTable()
    end if

    ! style buttons: show/hide
    if (iw_button("Show All##showallatoms")) then
       w%rep%atom_style(1:w%rep%natom_style)%shown = .true.
       changed = .true.
    end if
    call iw_tooltip("Show all atoms/bonds/labels in the system",ttshown)
    if (iw_button("Hide All##hideallatoms",sameline=.true.)) then
       w%rep%atom_style(1:w%rep%natom_style)%shown = .false.
       changed = .true.
    end if
    call iw_tooltip("Hide all atoms/bonds/labels in the system",ttshown)
    if (iw_button("Toggle Show/Hide##toggleallatoms",sameline=.true.)) then
       do i = 1, w%rep%natom_style
          w%rep%atom_style(i)%shown = .not.w%rep%atom_style(i)%shown
       end do
       changed = .true.
    end if
    call iw_tooltip("Toggle the show/hide status for all atoms/bonds/labels",ttshown)

    !!! atoms !!!

    str2 = "Atoms##atomsglobaldisplay " // c_null_char
    call igPushStyleColor_Vec4(ImGuiCol_Text,ColorHighlightText)
    changed = changed .or. igCheckbox(c_loc(str2),w%rep%atoms_display)
    call igPopStyleColor(1)
    call iw_tooltip("Draw the atoms in the scene",ttshown)

    if (w%rep%atoms_display) then
       ! style buttons: set radii
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
          do i = 1, w%rep%natom_style
             if (w%rep%atom_style_type == 0) then
                ! species
                ispc = i
             elseif (w%rep%atom_style_type == 1) then
                ! nneq
                ispc = sys(isys)%c%at(i)%is
             elseif (w%rep%atom_style_type == 2) then
                ! ncel
                ispc = sys(isys)%c%atcel(i)%is
             end if
             iz = sys(isys)%c%spc(ispc)%z
             if (w%rep%atom_radii_reset_type == 0) then
                w%rep%atom_style(i)%rad = real(atmcov(iz),c_float)
             else
                w%rep%atom_style(i)%rad = real(atmvdw(iz),c_float)
             end if
             w%rep%atom_style(i)%rad = w%rep%atom_style(i)%rad * w%rep%atom_radii_reset_scale
          end do
          changed = .true.
       end if

       ! style buttons: set color
       call iw_combo_simple("Colors ##atomcolorselect","jmol (light)" // c_null_char // "jmol2 (dark)" // c_null_char,&
          w%rep%atom_color_reset_type,changed=ch)
       call iw_tooltip("Set the color of all atoms to the tabulated values",ttshown)
       if (ch) then
          do i = 1, w%rep%natom_style
             if (w%rep%atom_style_type == 0) then
                ! species
                ispc = i
             elseif (w%rep%atom_style_type == 1) then
                ! nneq
                ispc = sys(isys)%c%at(i)%is
             elseif (w%rep%atom_style_type == 2) then
                ! ncel
                ispc = sys(isys)%c%atcel(i)%is
             end if
             iz = sys(isys)%c%spc(ispc)%z
             if (w%rep%atom_color_reset_type == 0) then
                w%rep%atom_style(i)%rgb = real(jmlcol(:,iz),c_float) / 255._c_float
             else
                w%rep%atom_style(i)%rgb = real(jmlcol2(:,iz),c_float) / 255._c_float
             end if
          end do
          changed = .true.
       end if
    end if

    !!! bonds !!!

    str2 = "Bonds##bondsglobaldisplay" // c_null_char
    call igPushStyleColor_Vec4(ImGuiCol_Text,ColorHighlightText)
    changed = changed .or. igCheckbox(c_loc(str2),w%rep%bonds_display)
    call igPopStyleColor(1)
    call iw_tooltip("Draw the bonds in the scene",ttshown)

    if (w%rep%bonds_display) then
       ! bond styles
       !! colors
       call iw_combo_simple("Style ##bondcolorsel","Single color"//c_null_char//"Two colors (atoms)"//c_null_char,&
          w%rep%bond_color_style,changed=ch)
       call iw_tooltip("Style for the bond colors",ttshown)
       if (ch) changed = .true.
       if (w%rep%bond_color_style == 0) then
          call igSameLine(0._c_float,-1._c_float)
          str2 = "##bondcolor" // c_null_char
          ch = igColorEdit3(c_loc(str2),w%rep%bond_rgb,ImGuiColorEditFlags_NoInputs)
          call iw_tooltip("Color for the representation bonds",ttshown)
          call iw_clamp_color3(w%rep%bond_rgb)
          call iw_text("Color",sameline=.true.)
          if (ch) then
             w%rep%bond_rgb = min(w%rep%bond_rgb,1._c_float)
             w%rep%bond_rgb = max(w%rep%bond_rgb,0._c_float)
             changed = .true.
          end if
       end if

       !! radius
       call igSameLine(0._c_float,-1._c_float)
       str2 = "Radius ##bondradius" // c_null_char
       str3 = "%.3f" // c_null_char
       call igPushItemWidth(iw_calcwidth(5,1))
       changed = changed .or. igDragFloat(c_loc(str2),w%rep%bond_rad,0.005_c_float,0._c_float,2._c_float,&
          c_loc(str3),ImGuiSliderFlags_AlwaysClamp)
       call igPopItemWidth()
       call iw_tooltip("Radii of the cylinders representing the bonds",ttshown)
    end if

    !!! labels !!!

    str2 = "Labels##labelsglobaldisplay" // c_null_char
    call igPushStyleColor_Vec4(ImGuiCol_Text,ColorHighlightText)
    changed = changed .or. igCheckbox(c_loc(str2),w%rep%labels_display)
    call igPopStyleColor(1)
    call iw_tooltip("Draw the atom labels in the scene",ttshown)

    if (w%rep%labels_display) then
       ! label styles
       if (sys(isys)%c%ismolecule) then
          lst = lsttrans(w%rep%label_style)
          call iw_combo_simple("Text##labelcontentselect","Atom name"// c_null_char// "Atom ID"// c_null_char//&
             "Species ID"// c_null_char// "Atomic number"// c_null_char// "Molecule ID"// c_null_char,&
             lst,changed=ch)
          w%rep%label_style = lsttransi(lst)
       else
          call iw_combo_simple("Text##labelcontentselect","Atom name"// c_null_char// "Cell atom ID"// c_null_char//&
             "Cell atom ID + lattice vector"// c_null_char// "Symmetry-unique atom ID"// c_null_char//&
             "Species ID"// c_null_char// "Atomic number"// c_null_char// "Molecule ID"// c_null_char,&
             w%rep%label_style,changed=ch)
       end if
       call iw_tooltip("Text to display in the atom labels",ttshown)
       changed = changed .or. ch

       ! scale and constant size
       str2 = "Scale##labelscale" // c_null_char
       str3 = "%.2f" // c_null_char
       call igPushItemWidth(iw_calcwidth(4,1))
       changed = changed .or. igDragFloat(c_loc(str2),w%rep%label_scale,0.01_c_float,0._c_float,10._c_float,c_loc(str3),&
          ImGuiSliderFlags_AlwaysClamp)
       call igPopItemWidth()
       call iw_tooltip("Scale factor for the atom labels",ttshown)

       call igSameLine(0._c_float,-1._c_float)
       str2 = "Constant size##labelconstsize" // c_null_char
       changed = changed .or. igCheckbox(c_loc(str2),w%rep%label_const_size)
       call iw_tooltip("Labels have constant size (on) or labels scale with the size of the associated atom (off)",ttshown)

       ! color
       call igSameLine(0._c_float,-1._c_float)
       str2 = "Color##labelcolor" // c_null_char
       changed = changed .or. igColorEdit3(c_loc(str2),w%rep%label_rgb,ImGuiColorEditFlags_NoInputs)
       call iw_tooltip("Color of the atom labels",ttshown)
       call iw_clamp_color3(w%rep%label_rgb)
    end if

  end function draw_editrep_atoms

  !> Draw the editrep window, unit cell class. Returns true if the
  !> scene needs rendering again. ttshown = the tooltip flag.
  module function draw_editrep_unitcell(w,ttshown) result(changed)
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_calcwidth, iw_radiobutton, iw_button,&
       iw_clamp_color3
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
    call iw_tooltip("Cell not repeated for this representation",ttshown)
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
    str1 = "Color crystallographic axes" // c_null_char
    changed = changed .or. igCheckbox(c_loc(str1),w%rep%uc_coloraxes)
    call iw_tooltip("Represent crystallographic axes with colors (a=red,b=green,c=blue)",ttshown)

    str1 = "Radius##outer" // c_null_char
    str2 = "%.3f" // c_null_char
    call igPushItemWidth(iw_calcwidth(5,1))
    changed = changed .or. igDragFloat(c_loc(str1),w%rep%uc_radius,0.005_c_float,0._c_float,&
       5._c_float,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
    w%rep%uc_radius = max(w%rep%uc_radius,0._c_float)
    call igPopItemWidth()
    call iw_tooltip("Radii of the unit cell edges",ttshown)

    call igSameLine(0._c_float,-1._c_float)
    str1 = "Color" // c_null_char
    ch = igColorEdit3(c_loc(str1),w%rep%uc_rgb,ImGuiColorEditFlags_NoInputs)
    call iw_tooltip("Color of the unit cell edges",ttshown)
    call iw_clamp_color3(w%rep%uc_rgb)
    if (ch) then
       w%rep%uc_rgb = min(w%rep%uc_rgb,1._c_float)
       w%rep%uc_rgb = max(w%rep%uc_rgb,0._c_float)
       changed = .true.
    end if

    !! inner divisions
    call iw_text("Inner Divisions",highlight=.true.)
    str1 = "Display inner divisions" // c_null_char
    changed = changed .or. igCheckbox(c_loc(str1),w%rep%uc_inner)
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

       str1 = "Use dashed lines" // c_null_char
       changed = changed .or. igCheckbox(c_loc(str1),w%rep%uc_innerstipple)
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
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG,&
       BIND_CLOSE_ALL_DIALOGS
    use tools_io, only: ferror, string
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
       w%okfile = ""
       w%okfilter = "   "
       w%nsample = 16
       w%jpgquality = 90
       w%exportview = .true.
       w%npixel = win(w%idparent)%FBOside
       w%transparentbg = .true.
       w%errmsg = ""
    end if

    ! initialize
    doquit = .false.
    isys = win(w%idparent)%view_selected

    ! check if we have info from the save image file window when it
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

    str2 = "Export the viewport only" // c_null_char
    ldum = igCheckbox(c_loc(str2),w%exportview)
    call iw_tooltip("Export the viewport only or the whole render buffer",ttshown)

    str2 = "Transparent background" // c_null_char
    ldum = igCheckbox(c_loc(str2),w%transparentbg)
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
       call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
       call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
       call glBindTexture(GL_TEXTURE_2D, 0)

       ! render buffer
       call glBindRenderbuffer(GL_RENDERBUFFER, msFBOdepth)
       call glRenderbufferStorageMultisample(GL_RENDERBUFFER, w%nsample, GL_DEPTH32F_STENCIL8, w%npixel, w%npixel)
       call glBindRenderbuffer(GL_RENDERBUFFER, 0)
       call glBindRenderbuffer(GL_RENDERBUFFER, endFBOdepth)
       call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH32F_STENCIL8, w%npixel, w%npixel)
       call glBindRenderbuffer(GL_RENDERBUFFER, 0)

       ! frame buffer
       call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
       call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, msFBOtex, 0)
       call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, msFBOdepth)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
          w%errmsg = "Failed creating the multi-sampled framebuffer (too large?)"
          goto 999
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, endFBO)
       call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, endFBOtex, 0)
       call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, endFBOdepth)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
          w%errmsg = "Failed creating the render framebuffer (too large?)"
          goto 999
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! render the scene to the multisampled framebuffer
       call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
       call glViewport(0_c_int,0_c_int,w%npixel,w%npixel)
       if (w%transparentbg) then
          call glClearColor(sysc(isys)%sc%bgcolor(1),sysc(isys)%sc%bgcolor(2),sysc(isys)%sc%bgcolor(3),0._c_float)
       else
          call glClearColor(sysc(isys)%sc%bgcolor(1),sysc(isys)%sc%bgcolor(2),sysc(isys)%sc%bgcolor(3),1._c_float)
       end if
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
       goodsys = (isys >= 1 .and. isys <= nsys)
       if (goodsys) goodsys = (sysc(isys)%status == sys_init)
       if (goodsys) call sysc(isys)%sc%render()
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! blit the multisampled buffer to the normal colorbuffer
       call glBindFramebuffer(GL_READ_FRAMEBUFFER, msFBO)
       call glBindFramebuffer(GL_DRAW_FRAMEBUFFER, endFBO)
       call glBlitFramebuffer(0, 0, w%npixel, w%npixel, 0, 0, w%npixel, w%npixel, GL_COLOR_BUFFER_BIT, GL_NEAREST)
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
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
          w%errmsg = "Could not save image to file: " // string(w%okfile)

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

end submodule view
