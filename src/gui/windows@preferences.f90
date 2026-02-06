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

! Routines for the preferences window.
submodule (windows) preferences
  use interfaces_cimgui
  implicit none
contains

  !> Draw the preferences window
  module subroutine draw_preferences(w)
    use windows, only: nwin, win
    use tools_io, only: nameguess, string
    use gui_main, only: g, tooltip_enabled, tooltip_delay, tooltip_wrap_factor,&
       tree_select_updates_inpcon, tree_select_updates_view, io,&
       set_default_ui_settings, ColorTableCellBg, ColorHighlightScene,&
       ColorHighlightSelectScene, ColorHighlightSelectScene, ColorMeasureSelect, &
       ColorElement
    use systems, only: nsys, sysc
    use interfaces_cimgui
    use keybindings
    use utils, only: iw_tooltip, iw_button, iw_text, iw_calcwidth, iw_clamp_color4,&
       iw_checkbox, iw_coloredit, iw_dragfloat_realc
    use param, only: maxzat0
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, str2, zeroc
    character(len=:), allocatable, target :: strf
    logical(c_bool) :: ldum, ch
    logical :: doquit
    type(ImVec2) :: sz, szero
    integer :: i, newkey, igroup, kmod, nrow, iwin, isys
    integer(c_int) :: flags
    real(c_float) :: width

    logical, save :: ttshown = .false. ! tooltip flag
    type(c_ptr), save :: cfilter = c_null_ptr ! filter object (allocated first pass, never destroyed)
    integer(c_int), save :: catid = 0 ! category ID (from left panel) 0=interface,1=keybinding,2=colors
    integer(c_int), save :: getbind = -1 ! get binding flag

    real(c_float), parameter :: wleft = 120._c_float

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! initialize
    doquit = .false.
    szero%x = 0
    szero%y = 0

    ! text filter
    call igAlignTextToFramePadding()
    call iw_text("Filter")
    call igSameLine(0._c_float,-1._c_float)
    if (.not.c_associated(cfilter)) &
       cfilter = ImGuiTextFilter_ImGuiTextFilter(c_loc(zeroc))
    str = "##preferencesfilter" // c_null_char
    ldum = ImGuiTextFilter_Draw(cfilter,c_loc(str),0._c_float)
    call iw_tooltip("Filter UI settings by name in the list below. Use comma-separated fields&
       & and - for excluding. Example: inc1,inc2,-exc includes all settings&
       & with inc1 or inc2 and excludes settings with exc.",ttshown)
    if (iw_button("Clear",sameline=.true.)) then
       if (c_associated(cfilter)) &
          call ImGuiTextFilter_Clear(cfilter)
    end if
    call iw_tooltip("Clear the filter",ttshown)

    ! left panel with selectables
    str = "leftpanel" // c_null_char
    sz%x = wleft
    sz%y = 0._c_float
    if (igBeginChild_Str(c_loc(str),sz,.true._c_bool,ImGuiWindowFlags_None)) then
       str = "Interface" // c_null_char
       if (igSelectable_Bool(c_loc(str),logical(catid == 0,c_bool),ImGuiSelectableFlags_None,szero)) catid = 0
       call iw_tooltip("Global settings for the graphical user interface",ttshown)
       str = "Key bindings" // c_null_char
       if (igSelectable_Bool(c_loc(str),logical(catid == 1,c_bool),ImGuiSelectableFlags_None,szero)) catid = 1
       call iw_tooltip("User interface key bindings",ttshown)
       str = "Colors" // c_null_char
       if (igSelectable_Bool(c_loc(str),logical(catid == 2,c_bool),ImGuiSelectableFlags_None,szero)) catid = 2
       call iw_tooltip("Background colors for the systems in the tree window",ttshown)
    end if
    call igEndChild()
    call igSameLine(0._c_float,-1._c_float)

    ! right panel
    call igBeginGroup()
    str = "rightpanel" // c_null_char
    sz%x = 0._c_float
    sz%y = -igGetFrameHeightWithSpacing() - g%Style%ItemSpacing%y
    if (igBeginChild_Str(c_loc(str),sz,.false._c_bool,ImGuiWindowFlags_None)) then

       if (catid == 0) then
          !! Interface
          call iw_text("Interface",highlight=.true.)
          call igSeparator()

          str = "Font scale" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_dragfloat_realc(str,x1=io%FontGlobalScale,speed=0.005_c_float,min=0.3_c_float,&
                max=3.0_c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Scale factor for the user interface font",ttshown)
          end if

          str = "Enable tooltips" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_checkbox(str,tooltip_enabled)
             call iw_tooltip("Show/hide the tooltips when hovering interface elements with the mouse",ttshown)
          end if

          str = "Tooltip delay (s)" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_dragfloat_realc(str,x1=tooltip_delay,speed=0.1_c_float,min=0._c_float,&
                max=5._c_float,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Delay for showing the tooltips",ttshown)
          end if

          str = "Tooltip maximum width (pixels)" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_dragfloat_realc(str,x1=tooltip_wrap_factor,speed=1._c_float,&
                min=0._c_float,max=1000._c_float,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
             call iw_tooltip("Width of the interface tooltips",ttshown)
          end if

          str = "Tree selects input console system" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_checkbox(str,tree_select_updates_inpcon)
             call iw_tooltip("Selecting a system on the tree changes the system in the input console",ttshown)
          end if

          str = "Tree selects view system" // c_null_char
          if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str),c_null_ptr)) then
             ldum = iw_checkbox(str,tree_select_updates_view)
             call iw_tooltip("Selecting a system on the tree changes the system in the view window",ttshown)
          end if

       elseif (catid == 1) then
          !! key bindings
          do igroup = 1, group_NUM
             call iw_text(trim(groupnames(igroup)),highlight=.true.)
             if (igroup == 1) then
                call iw_text("(?)",sameline=.true.)
                call iw_tooltip("Left click to assign a new binding. Right-click to toggle double click behavior&
                   & (only for mouse input). Middle click to erase the binding.")
             end if
             call igSeparator()

             str = "##keybindtable" // c_null_char
             flags = ImGuiTableFlags_NoSavedSettings
             flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
             flags = ior(flags,ImGuiTableFlags_NoBordersInBody)
             sz%x = 0
             sz%y = 0
             if (igBeginTable(c_loc(str),2,flags,sz,0._c_float)) then
                ! set up the columns
                width = iw_calcwidth(len(bindnames(1)),0)
                call igTableSetupColumn(c_null_ptr,ImGuiTableColumnFlags_None,width,0)
                call igTableSetupColumn(c_null_ptr,ImGuiTableColumnFlags_WidthStretch,0.0_c_float,1)
                call igTableSetColumnWidthAutoAll(igGetCurrentTable())

                ! table rows
                do i = 1, BIND_NUM

                   if (groupbind(i) /= igroup) cycle
                   str2 = trim(bindnames(i)) // c_null_char
                   if (.not.ImGuiTextFilter_PassFilter(cfilter,c_loc(str2),c_null_ptr)) cycle

                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   if (igTableSetColumnIndex(0)) then
                      call igAlignTextToFramePadding()
                      call iw_text(str2)
                   end if

                   if (igTableSetColumnIndex(1)) then
                      strf = trim(get_bind_keyname(i))

                      call igPushID_Int(int(i,c_int))
                      if (iw_button(strf)) getbind = i
                      call igPopID()

                      if (igIsItemHovered(ImGuiHoveredFlags_None)) then
                         ! right click to toggle
                         if (igIsMouseClicked(ImGuiPopupFlags_MouseButtonRight,.false._c_bool)) then
                            newkey = ImGuiKey_None
                            if (keybind(i) == ImGuiKey_MouseLeft) then
                               newkey = ImGuiKey_MouseLeftDouble
                            elseif (keybind(i) == ImGuiKey_MouseLeftDouble) then
                               newkey = ImGuiKey_MouseLeft
                            elseif (keybind(i) == ImGuiKey_MouseRight) then
                               newkey = ImGuiKey_MouseRightDouble
                            elseif (keybind(i) == ImGuiKey_MouseRightDouble) then
                               newkey = ImGuiKey_MouseRight
                            elseif (keybind(i) == ImGuiKey_MouseMiddle) then
                               newkey = ImGuiKey_MouseMiddleDouble
                            elseif (keybind(i) == ImGuiKey_MouseMiddleDouble) then
                               newkey = ImGuiKey_MouseMiddle
                            end if
                            if (newkey /= ImGuiKey_None) then
                               call set_bind(i,newkey,modbind(i))
                            end if
                         end if

                         ! middle click to erase
                         if (igIsMouseClicked(ImGuiPopupFlags_MouseButtonMiddle,.false._c_bool)) &
                            call set_bind(i,ImGuiKey_None,modbind(i))
                      end if
                   end if
                end do
                call igEndTable()
             end if
          end do

          ! popup to read a key
          if (getbind /= -1) then
             use_keybindings = .false.
             str2 = "Choosekey" // c_null_char
             flags = ImGuiWindowFlags_AlwaysAutoResize
             flags = ior(flags,ImGuiWindowFlags_NoTitleBar)
             flags = ior(flags,ImGuiWindowFlags_NoMove)
             call igOpenPopup_Str(c_loc(str2),ImGuiPopupFlags_None)
             if (igBeginPopupModal(c_loc(str2),logical(.true.,c_bool),flags)) then
                call iw_text("Please press a key or mouse button.")
                if (set_bind_from_user_input(getbind)) then
                   getbind = -1
                   use_keybindings = .true.
                   call igCloseCurrentPopup()
                end if
                call igEndPopup()
             end if
          end if

       elseif (catid == 2) then
          !! tree colors
          call iw_text("Tree window",highlight=.true.)
          call igSeparator()

          ! color editors
          call color_edit4("Crystal","A molecular system with a single molecule",ColorTableCellBg(:,0))
          call color_edit4("Layered crystal","A layered crystal, with two-dimensional periodic bonding networks",&
             ColorTableCellBg(:,1))
          call color_edit4("Crystal with 1D bonded chains",&
             "A crystal, with a one-dimensional periodic bonding networks",&
             ColorTableCellBg(:,2))
          call color_edit4("Molecular crystal","A molecular crystal",ColorTableCellBg(:,3))
          call color_edit4("Surface","A 2D surface",ColorTableCellBg(:,4))
          call color_edit4("Chain","An one-dimensional atomic structure (e.g. a polymer)",ColorTableCellBg(:,5))
          call color_edit4("Molecule in a box","A molecule in a large periodic box",ColorTableCellBg(:,6))
          call color_edit4("Molecule","A molecule",ColorTableCellBg(:,7))
          call color_edit4("Molecular cluster","A molecular cluster",ColorTableCellBg(:,8))

          ! atom selection and highlights
          call iw_text("Atom selection and highlights",highlight=.true.)
          call igSeparator()
          call color_edit4("Hovered atom","Atoms hovered by mouse in a table are&
             & shown in this color in the view window",ColorHighlightScene)
          call color_edit4("Selected atom","Atoms selected in the view or edit&
             & geometry windows",ColorHighlightSelectScene)
          call color_edit4("Measure (1st atom)","Atoms selected by double-click&
             & when measuring distances and angles",ColorMeasureSelect(:,1))
          call color_edit4("Measure (2nd atom)","Atoms selected by double-click&
             & when measuring distances and angles",ColorMeasureSelect(:,2))
          call color_edit4("Measure (3rd atom)","Atoms selected by double-click&
             & when measuring distances and angles",ColorMeasureSelect(:,3))
          call color_edit4("Measure (4th atom)","Atoms selected by double-click&
             & when measuring distances and angles",ColorMeasureSelect(:,4))

          ! elements
          call igAlignTextToFramePadding()
          call iw_text("Elements",highlight=.true.)
          ldum = iw_checkbox("Apply changes immediately",w%color_preferences_reset_reps,&
             sameline=.true.)
          call iw_tooltip("If checked, changes to the element colors are immediately&
             & applied to all atoms in all open systems.")
          call igSeparator()

          str = "##elementcolortable" // c_null_char
          flags = ImGuiTableFlags_NoSavedSettings
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_Borders)
          if (igBeginTable(c_loc(str),11,flags,szero,0._c_float)) then
             ! header
             call igTableSetupColumn(c_null_ptr,ImGuiTableColumnFlags_None,iw_calcwidth(2,0),0)
             width = iw_calcwidth(5,0)
             do i = 0, 9
                str = string(i) // c_null_char
                call igTableSetupColumn(c_loc(str),ImGuiTableColumnFlags_None,width,i+1)
             end do
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())
             call igTableHeadersRow()

             nrow = -1
             do i = 0, maxzat0
                kmod = mod(i,10)
                if (kmod == 0) then
                   call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                   if (igTableSetColumnIndex(0)) then
                      call igAlignTextToFramePadding()
                      nrow = nrow + 1
                      call iw_text(string(nrow))
                   end if
                end if
                if (igTableSetColumnIndex(kmod+1)) then
                   ch = iw_coloredit(nameguess(i,.true.),rgb=ColorElement(:,i))
                   if (ch .and. w%color_preferences_reset_reps) then
                      !! reset colors
                      ! systems
                      do isys = 1, nsys
                         call sysc(isys)%sc%reset_atom_colors()
                      end do
                      ! alternate view windows
                      do iwin = 1, nwin
                         if (win(iwin)%type == wintype_view.and..not.win(iwin)%ismain.and.&
                            associated(win(iwin)%sc)) then
                            call win(iwin)%sc%reset_atom_colors()
                         end if
                      end do
                   end if
                end if
             end do

             call igEndTable()
          end if

       end if
    end if
    call igEndChild()

    ! child for final buttons
    str = "finalbuttons" // c_null_char
    sz%x = 0._c_float
    sz%y = 0._c_float
    if (igBeginChild_Str(c_loc(str),sz,.false._c_bool,ImGuiWindowFlags_None)) then
       call igSetCursorPosX(iw_calcwidth(10,2,from_end=.true.) - g%Style%ScrollbarSize)
       if (iw_button("Reset",danger=.true.)) &
          call set_default_ui_settings()
       call iw_tooltip("Reset to the default settings",ttshown)
       if (iw_button("Close",sameline=.true.)) doquit = .true.
       call iw_tooltip("Close this window",ttshown)
    end if
    call igEndChild()
    call igEndGroup()

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
       if (c_associated(cfilter)) call ImGuiTextFilter_Clear(cfilter)
       catid = 0_c_int
       getbind = -1_c_int
    end subroutine init_state

    ! terminate the state for this window
    subroutine end_state()
       if (c_associated(cfilter)) call ImGuiTextFilter_Clear(cfilter)
    end subroutine end_state

    ! display the color editor
    subroutine color_edit4(str,strtt,rgba)
      character*(*), intent(in) :: str, strtt
      real(c_float), intent(inout) :: rgba(4)

      character(len=:), allocatable, target :: str_
      logical(c_bool) :: ldum

      str_ = str // c_null_char
      if (ImGuiTextFilter_PassFilter(cfilter,c_loc(str_),c_null_ptr)) then
         ldum = igColorEdit4(c_loc(str_),rgba,ImGuiColorEditFlags_NoInputs)
         call iw_tooltip(strtt,ttshown)
         call iw_clamp_color4(rgba)
      end if

    end subroutine color_edit4

  end subroutine draw_preferences

end submodule preferences
