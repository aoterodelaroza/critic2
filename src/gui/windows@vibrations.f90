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

! Routines for the vibrations window.
submodule (windows) vibrations
  use interfaces_cimgui
  implicit none
contains

  !> Draw the vibrations window
  module subroutine draw_vibrations(w)
    use interfaces_glfw, only: glfwGetTime
    use crystalseedmod, only: crystalseed
    use scenes, only: anim_speed_default, anim_amplitude_default, anim_amplitude_max,&
       anim_speed_max
    use systems, only: sysc, sys, sys_init, add_systems_from_seeds,&
       launch_initialization_thread, ok_system
    use gui_main, only: g
    use utils, only: iw_text, iw_button, iw_tooltip, iw_calcheight, iw_calcwidth,&
       iw_combo_simple, iw_radiobutton, iw_dragfloat_realc
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use tools_math, only: rational_approx
    use tools_io, only: string, ioj_right
    use param, only: cm1tothz, bohrtoa
    class(window), intent(inout), target :: w

    logical(c_bool) :: selected
    logical :: doquit, goodsys, vib_ok, goodparent, ldum, fset
    integer :: isys, i, digits, iaux
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

    ! vibrations ok?
    goodsys = ok_system(isys,sys_init)
    vib_ok = goodsys
    if (vib_ok) vib_ok = sys(isys)%c%vib%hasvibs
    if (vib_ok) vib_ok = (sys(isys)%c%vib%nqpt > 0) .and. (sys(isys)%c%vib%nfreq > 0)
    if (vib_ok) vib_ok = associated(win(w%idparent)%sc)

    ! header
    if (goodsys) then
       ! system name
       call iw_text("System",highlight=.true.)
       call iw_text("(" // string(isys) // ") " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! source of vibration data
       call igAlignTextToFramePadding()
       call iw_text("Vibration data",highlight=.true.)

       if (iw_button("Clear",sameline=.true.,danger=.true.)) then
          call sys(isys)%c%vib%end()
          win(w%idparent)%sc%iqpt_selected = 0
          win(w%idparent)%sc%ifreq_selected = 0
          vib_ok = .false.
       end if
       call iw_tooltip("Clear vibration data for this system",ttshown)

       if (.not.vib_ok) then
          if (iw_button("Load",danger=.true.,sameline=.true.)) &
             iaux = stack_create_window(wintype_dialog,.true.,purpose=wpurp_dialog_openvibfile,isys=isys,orraise=-1)
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
          flags = ior(flags,ImGuiTableFlags_RowBg)
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
       flags = ior(flags,ImGuiTableFlags_RowBg)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       flags = ior(flags,ImGuiTableFlags_ScrollX)
       str1="##tablevibrationfreqs" // c_null_char
       sz0%x = 0._c_float
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
             call igPushItemWidth(iw_calcwidth(5,1))
             if (iw_dragfloat_realc("Amplitude##amplitude",x1=win(w%idparent)%sc%anim_amplitude,&
                speed=0.01_c_float,min=0._c_float,max=anim_amplitude_max,sformat="%.2f",flags=ImGuiSliderFlags_AlwaysClamp))&
                win(w%idparent)%forcerender = .true.
             call igPopItemWidth()
             call iw_tooltip("Amplitude of the atomic displacements",ttshown)

             call igSameLine(0._c_float,-1._c_float)
             call igPushItemWidth(iw_calcwidth(6,1))
             if (iw_dragfloat_realc("Phase##phase",x1=win(w%idparent)%sc%anim_phase,speed=0.001_c_float,&
                min=-1._c_float,max=1._c_float,sformat="%.3f",flags=ImGuiSliderFlags_AlwaysClamp)) &
                win(w%idparent)%forcerender = .true.
             call igPopItemWidth()
             call iw_tooltip("Phase for the atomic displacements along the chosen phonon normal mode",ttshown)
          else
             ! molecules
             call igPushItemWidth(iw_calcwidth(5,1))
             if (iw_dragfloat_realc("Displacement##amplitude",x1=win(w%idparent)%sc%anim_amplitude,&
                speed=0.01_c_float,min=-anim_amplitude_max,max=anim_amplitude_max,sformat="%.2f",&
                flags=ImGuiSliderFlags_AlwaysClamp))&
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
          call igPushItemWidth(iw_calcwidth(5,1))
          ldum = iw_dragfloat_realc("Amplitude##amplitude",x1=win(w%idparent)%sc%anim_amplitude,&
             speed=0.01_c_float,min=0._c_float,max=anim_amplitude_max,sformat="%.2f",flags=ImGuiSliderFlags_AlwaysClamp)
          call igPopItemWidth()
          call iw_tooltip("Amplitude of the atomic displacements",ttshown)

          call igSameLine(0._c_float,-1._c_float)
          call igPushItemWidth(iw_calcwidth(5,1))
          if (iw_dragfloat_realc("Speed##speed",x1=win(w%idparent)%sc%anim_speed,&
             speed=0.02_c_float,min=0.0_c_float,max=anim_speed_max,sformat="%.2f",flags=ImGuiSliderFlags_AlwaysClamp)) &
             win(w%idparent)%sc%timerefanimation = glfwGetTime()
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
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit = close the window
    if (doquit) &
       call w%end()

  end subroutine draw_vibrations

end submodule vibrations
