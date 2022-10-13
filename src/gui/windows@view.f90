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

  !> Draw the view.
  module subroutine draw_view(w)
    use interfaces_opengl3
    use interfaces_cimgui
    use gui_main, only: sysc, sys_init, nsys
    use utils, only: iw_text, iw_button, iw_tooltip
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer :: i
    type(ImVec2) :: szavail, sz0, sz1, szero
    type(ImVec4) :: tint_col, border_col
    character(kind=c_char,len=:), allocatable, target :: str1, str2
    logical(c_bool) :: ldum, is_selected
    logical :: hover
    integer(c_int) :: amax
    real(c_float) :: scal

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize
    szero%x = 0
    szero%y = 0

    ! gear menu
    str1="##viewgear" // c_null_char
    if (iw_button("⚙")) then
       call igOpenPopup_Str(c_loc(str1),ImGuiPopupFlags_None)
    end if
    if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_None)) then
       str2 = "Empty for now" // c_null_char
       ldum = igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)
       call igEndPopup()
    end if
    call iw_tooltip("Change the view options")

    ! the selected system combo
    call igSameLine(0._c_float,-1._c_float)
    str2 = "" // c_null_char
    if (w%view_selected >= 1 .and. w%view_selected <= nsys) then
       if (sysc(w%view_selected)%status == sys_init) &
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
             end if
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Choose the system displayed in the view",ttshown)

    ! get the remaining size for the texture
    call igGetContentRegionAvail(szavail)

    ! resize the render texture if not large enough
    amax = max(ceiling(max(szavail%x,szavail%y)),1)
    if (amax > w%FBOside) then
       amax = max(ceiling(1.5 * ceiling(max(szavail%x,szavail%y))),1)
       call w%delete_texture_view()
       call w%create_texture_view(amax)
       w%forcerender = .true.
    end if

    ! render the image to the texture, if requested
    if (w%forcerender) then
       call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
       call glViewport(0_c_int,0_c_int,w%FBOside,w%FBOside)
       call glClearColor(0.,0.,0.,0.)
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))

       if (w%view_selected > 0 .and. w%view_selected <= nsys) &
          call sysc(w%view_selected)%sc%render()

       call glBindFramebuffer(GL_FRAMEBUFFER, 0)
       w%forcerender = .false.
    end if

    ! draw the texture, largest region with the same shape as the available region
    ! that fits into the texture square
    scal = real(w%FBOside,c_float) / max(max(szavail%x,szavail%y),1._c_float)
    sz0%x = 0.5 * (real(w%FBOside,c_float) - szavail%x * scal) / real(w%FBOside,c_float)
    sz0%y = 0.5 * (real(w%FBOside,c_float) - szavail%y * scal) / real(w%FBOside,c_float)
    sz1%x = 1._c_float - sz0%x
    sz1%y = 1._c_float - sz0%y

    ! border and tint for the image, draw the image, update the rectangle
    tint_col%x = 1._c_float
    tint_col%y = 1._c_float
    tint_col%z = 1._c_float
    tint_col%w = 1._c_float
    border_col%x = 0._c_float
    border_col%y = 0._c_float
    border_col%z = 0._c_float
    border_col%w = 0._c_float
    call igImage(w%FBOtex, szavail, sz0, sz1, tint_col, border_col)
    hover = igIsItemHovered(ImGuiHoveredFlags_None)
    call igGetItemRectMin(w%v_rmin)
    call igGetItemRectMax(w%v_rmax)

    ! Process mouse events
    call w%process_events_view(hover)

  end subroutine draw_view

  !> Create the texture for the view window, with atex x atex pixels.
  module subroutine create_texture_view(w,atex)
    use interfaces_opengl3
    use tools_io, only: ferror, faterr
    class(window), intent(inout), target :: w
    integer, intent(in) :: atex

    ! FBO and buffers
    call glGenTextures(1, c_loc(w%FBOtex))
    call glGenRenderbuffers(1, c_loc(w%FBOdepth))
    call glGenFramebuffers(1, c_loc(w%FBO))

    ! texture
    call glBindTexture(GL_TEXTURE_2D, w%FBOtex)
    call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, atex, atex,&
       0, GL_RGBA, GL_UNSIGNED_BYTE, c_null_ptr)
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
    call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
    call glBindTexture(GL_TEXTURE_2D, 0)

    ! render buffer
    call glBindRenderbuffer(GL_RENDERBUFFER, w%FBOdepth)
    call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, atex, atex)
    call glBindRenderbuffer(GL_RENDERBUFFER, 0)

    ! frame buffer
    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
    call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, w%FBOtex, 0)
    call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, w%FBOdepth)

    ! check for errors
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
       call ferror('window_init','framebuffer is not complete',faterr)

    ! finish and write the texture size
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)
    w%FBOside = atex

  end subroutine create_texture_view

  !> Delete the texture for the view window
  module subroutine delete_texture_view(w)
    use interfaces_opengl3
    class(window), intent(inout), target :: w

    call glDeleteTextures(1, c_loc(w%FBOtex))
    call glDeleteRenderbuffers(1, c_loc(w%FBOdepth))
    call glDeleteFramebuffers(1, c_loc(w%FBO))

  end subroutine delete_texture_view

  !> Process the mouse events in the view window
  module subroutine process_events_view(w,hover)
    use interfaces_cimgui
    use scenes, only: scene
    use utils, only: translate, rotate
    use tools_math, only: cross_cfloat, matinv_cfloat
    use keybindings, only: is_bind_event, is_bind_mousescroll, BIND_NAV_ROTATE,&
       BIND_NAV_TRANSLATE, BIND_NAV_ZOOM, BIND_NAV_RESET
    use gui_main, only: io, nsys, sysc
    class(window), intent(inout), target :: w
    logical, intent(in) :: hover

    type(ImVec2) :: texpos, mousepos
    real(c_float) :: ratio, depth, pos3(3), vnew(3), vold(3), axis(3), lax
    real(c_float) :: mpos2(2), ang
    type(scene), pointer :: sc

    integer, parameter :: ilock_no = 1
    integer, parameter :: ilock_left = 2
    integer, parameter :: ilock_mid = 3
    integer, parameter :: ilock_right = 4
    integer, parameter :: ilock_scroll = 5
    integer, parameter :: ilock_NUM = 5

    real(c_float), parameter :: mousesens_zoom0 = 0.15_c_float
    real(c_float), parameter :: mousesens_rot0 = 3._c_float
    real(c_float), parameter :: min_zoom = 1._c_float
    real(c_float), parameter :: max_zoom = 100._c_float

    type(ImVec2), save :: mposlast
    real(c_float), save :: mpos0_r(3), mpos0_l(3), cpos0_r(3), cpos0_l(3), world0(4,4)
    real(c_float), save :: world0inv(3,3), mpos0_s
    integer, save :: ilock = ilock_no

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! only process if there is an associated system is viewed and scene is initialized
    if (w%view_selected < 1 .or. w%view_selected > nsys) return
    sc => sysc(w%view_selected)%sc
    if (.not.sc%isinit) return

    ! process global events
    ! if (hover && IsBindEvent(BIND_VIEW_ALIGN_A_AXIS,false))
    !   outb |= sc->alignViewAxis(1);
    ! else if (hover && IsBindEvent(BIND_VIEW_ALIGN_B_AXIS,false))
    !   outb |= sc->alignViewAxis(2);
    ! else if (hover && IsBindEvent(BIND_VIEW_ALIGN_C_AXIS,false))
    !   outb |= sc->alignViewAxis(3);
    ! else if (hover && IsBindEvent(BIND_VIEW_ALIGN_X_AXIS,false))
    !   outb |= sc->alignViewAxis(-1);
    ! else if (hover && IsBindEvent(BIND_VIEW_ALIGN_Y_AXIS,false))
    !   outb |= sc->alignViewAxis(-2);
    ! else if (hover && IsBindEvent(BIND_VIEW_ALIGN_Z_AXIS,false))
    !   outb |= sc->alignViewAxis(-3);

    ! process mode-specific events
    if (w%view_mousebehavior == MB_Navigation) then
       call igGetMousePos(mousepos)
       texpos = mousepos

       ! transform to the texture pos
       if (abs(texpos%x) < 1e20 .and. abs(texpos%y) < 1e20) then
          call w%mousepos_to_texpos(texpos)
       end if
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

          pos3 = sc%campos - sc%scenecenter
          pos3 = pos3 - ratio * pos3
          if (norm2(pos3) < min_zoom) &
             pos3 = pos3 / norm2(pos3) * min_zoom
          if (norm2(pos3) > max_zoom * sc%scenerad) &
             pos3 = pos3 / norm2(pos3) * (max_zoom * sc%scenerad)
          sc%campos = sc%scenecenter + pos3

          call sc%update_view_matrix()
          call sc%update_projection_matrix()
          w%forcerender = .true.
       end if

       ! drag
       if (hover .and. is_bind_event(BIND_NAV_TRANSLATE,.false.) .and. (ilock == ilock_no .or. ilock == ilock_right)) then
          depth = w%texpos_viewdepth(texpos)
          if (depth < 1._c_float) then
             mpos0_r = (/texpos%x,texpos%y,depth/)
          else
             pos3 = 0._c_float
             call w%view_to_texpos(pos3)
             mpos0_r = (/texpos%x,texpos%y,pos3(3)/)
          end if
          cpos0_r = (/sc%campos(1),sc%campos(2),zero/)

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

                sc%campos(1) = cpos0_r(1) - (vnew(1) - vold(1))
                sc%campos(2) = cpos0_r(2) - (vnew(2) - vold(2))
                call sc%update_view_matrix()

                mposlast = mousepos
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
                   world0inv = sc%world(1:3,1:3)
                   call matinv_cfloat(world0inv,3)
                   axis = matmul(world0inv,axis)
                   mpos2(1) = texpos%x - mpos0_l(1)
                   mpos2(2) = texpos%y - mpos0_l(2)
                   ang = 2._c_float * norm2(mpos2) * mousesens_rot0 / w%FBOside

                   sc%world = translate(sc%world,real(sc%scenecenter,c_float))
                   sc%world = rotate(sc%world,ang,axis)
                   sc%world = translate(sc%world,real(-sc%scenecenter,c_float))

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
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      mposlast%x = 0._c_float
      mposlast%y = 0._c_float
      mpos0_r = 0._c_float
      mpos0_l = 0._c_float
      cpos0_r = 0._c_float
      cpos0_l = 0._c_float
      world0 = 0._c_float
      world0inv = 0._c_float
      mpos0_s = 0._c_float
      ilock = ilock_no
    end subroutine init_state
    subroutine end_state()
    end subroutine end_state

  end subroutine process_events_view

  !> Mouse position to texture position (screen coordinates)
  module subroutine mousepos_to_texpos(w,pos)
    class(window), intent(inout), target :: w
    type(ImVec2), intent(inout) :: pos

    real(c_float) :: dx, dy, xratio, yratio

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
  module function texpos_viewdepth(w,pos)
    use interfaces_opengl3
    class(window), intent(inout), target :: w
    type(ImVec2), intent(inout) :: pos
    real(c_float) :: texpos_viewdepth

    real(c_float), target :: depth

    ! pixels have the origin (0,0) at top left; the other end is bottom right, (FBOside-1,FBOside-1)
    call glBindFramebuffer(GL_FRAMEBUFFER, w%FBO)
    call glReadPixels(int(pos%x), int(w%FBOside - pos%y), 1_c_int, 1_c_int, GL_DEPTH_COMPONENT, GL_FLOAT, c_loc(depth))
    call glBindFramebuffer(GL_FRAMEBUFFER, 0)
    texpos_viewdepth = depth

  end function texpos_viewdepth

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

end submodule view
